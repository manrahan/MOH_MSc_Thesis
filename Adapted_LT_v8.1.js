/**
 * @license
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 * @author Justin Braaten (Google)
 * @author Zhiqiang Yang (USDA Forest Service)
 * @author Robert Kennedy (Oregon State University)
 * MODIFIED 1/2022 Ben Roberts-Pierel (Oregon State University)
 * 
 * @description This file contains functions for working with the LandTrendr
 * change detection algorithm in Google Earth Engine. For information on
 * LandTrendr and usage of functions in this file see
 * https://github.com/eMapR/LT-GEE. Please post issues to
 * https://github.com/eMapR/LT-GEE/issues.
 */

// #############################################################################
// ### VERSION ###
// #############################################################################

exports.version = 'LT: 0.2.0 ... MOH adaptation: 0.2';


print('IMPORTANT! Please be advised:');
print('- This version of the Adapted_LT.js modules') 
print('  uses some code adapted from the aut/or: @author Justin Braaten (Google) * @author Zhiqiang Yang (USDA Forest Service) * @author Robert Kennedy (Oregon State University)')
print('The latest edits to this code occur: 08/03/2023 for the adaptation efforts by @Mike OHanrahan (TU DELFT MSc research)')

//########################################################################################################
//##### ANNUAL SR TIME SERIES COLLECTION BUILDING FUNCTIONS ##### 
//########################################################################################################

//------ FILTER A COLLECTION FUNCTION -----
var filterCollection = function(year, startDay, endDay, sensor, aoi){
  return ee.ImageCollection('LANDSAT/'+ sensor + '/C02/T1_L2')
           .filterBounds(aoi)
           .filterDate(year+'-'+startDay, year+'-'+endDay);
};


//------ BUILD A COLLECTION FOR A GIVEN SENSOR AND YEAR -----
var buildSensorYearCollection = function(year, startDay, endDay, sensor, aoi, exclude){
  var startMonth = parseInt(startDay.substring(0, 2));
  var endMonth = parseInt(endDay.substring(0, 2));
  var srCollection;
  if(startMonth > endMonth){
    var oldYear = (parseInt(year)-1).toString();
    var newYear = year;
    var oldYearStartDay = startDay;
    var oldYearEndDay = '12-31';
    var newYearStartDay = '01-01';
    var newYearEndDay = endDay;
    
    var oldYearCollection = filterCollection(oldYear, oldYearStartDay, oldYearEndDay, sensor, aoi);
    var newYearCollection = filterCollection(newYear, newYearStartDay, newYearEndDay, sensor, aoi);
    
    srCollection = ee.ImageCollection(oldYearCollection.merge(newYearCollection));
  } else {
    srCollection = filterCollection(year, startDay, endDay, sensor, aoi);
  }
  
  srCollection = removeImages(srCollection, exclude)
  
  
  return srCollection;
};
exports.buildSensorYearCollection = buildSensorYearCollection


//------ RETRIEVE A SENSOR SR COLLECTION FUNCTION -----
//scaling values source: https://www.usgs.gov/faqs/how-do-i-use-scale-factor-landsat-level-2-science-products
//define a function to apply Collection 2 scaling coefficients 

var scaleLTdata = function(img){ 
  return ((img.multiply(0.0000275)).add(-0.2)).multiply(10000).toUint16();
}; 

var getSRcollection = function(year, startDay, endDay, sensor, aoi, maskThese, exclude) {
  // make sure that mask labels are correct
  maskThese = (typeof maskThese !== 'undefined') ?  maskThese : ['cloud','shadow','snow','water'];
  //var maskOptions = ['cloud', 'shadow', 'snow', 'water'];
  var maskOptions = ['cloud', 'shadow', 'snow', 'water', 'waterplus','nonforest']; // add new water and forest mask here Peter Clary 5/20/2020
  for(var i in maskThese){
    maskThese[i] =  maskThese[i].toLowerCase();
    var test = maskOptions.indexOf(maskThese[i]);
    if(test == -1){
      print('error: '+maskThese[i]+' is not included in the list maskable features. Please see ___ for list of maskable features to include in the maskThese parameter');
      return 'error';
    }
  }
  
  // get a landsat collection for given year, day range, and sensor
  var srCollection = buildSensorYearCollection(year, startDay, endDay, sensor, aoi, exclude);
  // apply the harmonization function to LC08 (if LC08), subset bands, unmask, and resample           
  srCollection = srCollection.map(function(img) {
    var dat = ee.Image(
      ee.Algorithms.If(
        (sensor == 'LC08') || (sensor == 'LC09'),                            // condition - if image is OLI
        scaleLTdata(img.select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'],['B1', 'B2', 'B3', 'B4', 'B5', 'B7'])).unmask(),

        scaleLTdata(img.select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7'],['B1', 'B2', 'B3', 'B4', 'B5', 'B7']))                   // false - else select out the reflectance bands from the non-OLI image
           .unmask()                                                       // ...unmask any previously masked pixels 
           //.resample('bicubic')                                          // ...resample by bicubic 
           .set('system:time_start', img.get('system:time_start'))         // ...set the output system:time_start metadata to the input image time_start otherwise it is null
      )
    );
    
    // makes a global forest mask
    var forCol = ee.ImageCollection("COPERNICUS/Landcover/100m/Proba-V/Global"); //PETER ADD
    var imgFor = forCol.toBands(); //PETER ADD
    var forestimage = imgFor.select('2015_forest_type') //PETER ADD
    
    // Computes the forest mask into a binary using an expression.
    var selectedForests = forestimage.expression( //PETER ADD
        'Band >= 0 ? 1 : 0', { //PETER ADD
          'Band': forestimage //PETER ADD
    }).clip(aoi); //PETER ADD
    
    //makes a global water mask
    var MappedWater = ee.Image("JRC/GSW1_1/GlobalSurfaceWater"); //PETER ADD
    // calculates water persistence 0 to 100 //PETER ADD
    var MappedWaterBinary = MappedWater.expression( //PETER ADD 
      'band > 99 ? 0 :  1  ', { //PETER ADD
        'band': MappedWater.select('recurrence') //PETER ADD
    }).clip(aoi); //PETER ADD
    
    
    var mask = ee.Image(1);
    if(maskThese.length !== 0){
      var qa = img.select('QA_PIXEL'); 
      for(var i in maskThese){
        if(maskThese[i] == 'water'){mask = qa.bitwiseAnd(1<<7).eq(0).multiply(mask)}
        if(maskThese[i] == 'shadow'){mask = qa.bitwiseAnd(1<<4).eq(0).multiply(mask)} 
        if(maskThese[i] == 'snow'){mask = qa.bitwiseAnd(1<<5).eq(0).multiply(mask)}
        if(maskThese[i] == 'cloud'){mask = qa.bitwiseAnd(1<<3).eq(0).multiply(mask)} 
        if(maskThese[i] == 'waterplus'){mask = mask.mask(MappedWaterBinary)} //PETER ADD
        if(maskThese[i] == 'nonforest'){mask = mask.mask(selectedForests)} // PETER ADD
      }
      return dat.mask(mask); //apply the mask - 0's in mask will be excluded from computation and set to opacity=0 in display
    } else{
      return dat;
    }
  });

  return srCollection; // return the prepared collection
};
exports.getSRcollection = getSRcollection;


//------ FUNCTION TO COMBINE LT05, LE07, LC08 and LC09 COLLECTIONS -----
var getCombinedSRcollectionOrig = function(year, startDay, endDay, aoi, maskThese) {
    var lt5 = getSRcollection(year, startDay, endDay, 'LT05', aoi, maskThese);       // get TM collection for a given year, date range, and area
    var le7 = getSRcollection(year, startDay, endDay, 'LE07', aoi, maskThese);       // get ETM+ collection for a given year, date range, and area
    var lc8 = getSRcollection(year, startDay, endDay, 'LC08', aoi, maskThese);       // get OLI collection for a given year, date range, and area
    var lc9 = getSRcollection(year, startDay, endDay, 'LC09', aoi, maskThese);       // get OLI collection for a given year, date range, and area
    var mergedCollection = ee.ImageCollection(lt5.merge(le7).merge(lc8).merge(lc9)); // merge the individual sensor collections into one imageCollection object
    return mergedCollection;                                              // return the Imagecollection
};
var getCombinedSRcollection = function(year, startDay, endDay, aoi, maskThese, exclude) {
    exclude = (typeof exclude !== 'undefined') ?  exclude : {}; // default to not exclude any images
    var lt5 = getSRcollection(year, startDay, endDay, 'LT05', aoi, maskThese, exclude);       // get TM collection for a given year, date range, and area
    var le7 = getSRcollection(year, startDay, endDay, 'LE07', aoi, maskThese, exclude);       // get ETM+ collection for a given year, date range, and area
    var lc8 = getSRcollection(year, startDay, endDay, 'LC08', aoi, maskThese, exclude);       // get OLI collection for a given year, date range, and area
    var lc9 = getSRcollection(year, startDay, endDay, 'LC09', aoi, maskThese, exclude);       // get OLI collection for a given year, date range, and area
    var mergedCollection = ee.ImageCollection(lt5.merge(le7).merge(lc8).merge(lc9)); // merge the individual sensor collections into one imageCollection object
    //mergedCollection = removeImages(mergedCollection, exclude);
    return mergedCollection;                                              // return the Imagecollection
};
exports.getCombinedSRcollection = getCombinedSRcollection; 



//------ FUNCTION TO REDUCE COLLECTION TO SINGLE IMAGE PER YEAR BY MEDOID -----
/*
  LT expects only a single image per year in a time series, there are lost of ways to
  do best available pixel compositing - we have found that a mediod composite requires little logic
  is robust, and fast
  
  Medoids are representative objects of a data set or a cluster with a data set whose average 
  dissimilarity to all the objects in the cluster is minimal. Medoids are similar in concept to 
  means or centroids, but medoids are always members of the data set.
*/

// make a medoid composite with equal weight among indices
var medoidMosaic = function(inCollection, dummyCollection) {
  
  // fill in missing years with the dummy collection
  var imageCount = inCollection.toList(1).length();                                                            // get the number of images 
  var finalCollection = ee.ImageCollection(ee.Algorithms.If(imageCount.gt(0), inCollection, dummyCollection)); // if the number of images in this year is 0, then use the dummy collection, otherwise use the SR collection
  
  // calculate median across images in collection per band
  var median = finalCollection.median();                                                                       // calculate the median of the annual image collection - returns a single 6 band image - the collection median per band
  
  // calculate the different between the median and the observation per image per band
  var difFromMedian = finalCollection.map(function(img) {
    var diff = ee.Image(img).subtract(median).pow(ee.Image.constant(2));                                       // get the difference between each image/band and the corresponding band median and take to power of 2 to make negatives positive and make greater differences weight more
    return diff.reduce('sum').addBands(img);                                                                   // per image in collection, sum the powered difference across the bands - set this as the first band add the SR bands to it - now a 7 band image collection
  });
  
  // get the medoid by selecting the image pixel with the smallest difference between median and observation per band 
  return ee.ImageCollection(difFromMedian).reduce(ee.Reducer.min(7)).select([1,2,3,4,5,6], ['B1','B2','B3','B4','B5','B7']); // find the powered difference that is the least - what image object is the closest to the median of teh collection - and then subset the SR bands and name them - leave behind the powered difference band
};


//------ FUNCTION TO APPLY MEDOID COMPOSITING FUNCTION TO A COLLECTION -------------------------------------------
var buildMosaic = function(year, startDay, endDay, aoi, dummyCollection, maskThese, exclude) {                                                                      // create a temp variable to hold the upcoming annual mosiac
  exclude = (typeof exclude !== 'undefined') ?  exclude : {}; // default to not exclude any images
  var collection = getCombinedSRcollection(year, startDay, endDay, aoi, maskThese, exclude);  // get the SR collection
  var img = medoidMosaic(collection, dummyCollection)                     // apply the medoidMosaic function to reduce the collection to single image per year by medoid 
              .set('system:time_start', (new Date(year,8,1)).valueOf());  // add the year to each medoid image - the data is hard-coded Aug 1st 
  return ee.Image(img).toUint16();                                                   // return as image object
};


//------ FUNCTION TO BUILD ANNUAL MOSAIC COLLECTION ------------------------------
var buildSRcollection = function(startYear, endYear, startDay, endDay, aoi, maskThese, exclude) {
  exclude = (typeof exclude !== 'undefined') ?  exclude : {}; // default to not exclude any images
  var dummyCollection = ee.ImageCollection([ee.Image([0,0,0,0,0,0]).mask(ee.Image(0))]); // make an image collection from an image with 6 bands all set to 0 and then make them masked values
  var imgs = [];                                                                         // create empty array to fill
  for (var i = startYear; i <= endYear; i++) {                                           // for each year from hard defined start to end build medoid composite and then add to empty img array
    var tmp = buildMosaic(i, startDay, endDay, aoi, dummyCollection, maskThese, exclude);                    // build the medoid mosaic for a given year
    imgs = imgs.concat(tmp.set('system:time_start', (new Date(i,8,1)).valueOf()));       // concatenate the annual image medoid to the collection (img) and set the date of the image - hard coded to the year that is being worked on for Aug 1st
  }
  return ee.ImageCollection(imgs);                                                       // return the array img array as an image collection
};
exports.buildSRcollection = buildSRcollection;

//------ FUNCTION TO RETURN A LIST OF IMAGES THAT GO INTO ANNUAL SR COMPOSITE COLLECTION ------------------------------
function getImgID(img){return ee.String(ee.Image(img).get('system:id'));}
function getImgIndex(img){return ee.String(ee.Image(img).get('system:index'));}


var getCollectionIDlist = function(startYear, endYear, startDay, endDay, aoi, exclude) {
  exclude = (typeof exclude !== 'undefined') ?  exclude : {}; // default to not exclude any images
  var first = true;
  for (var i = startYear; i <= endYear; i++){
    var lt5 = buildSensorYearCollection(i, startDay, endDay, 'LT05', aoi, exclude);
    var le7 = buildSensorYearCollection(i, startDay, endDay, 'LE07', aoi, exclude);
    var lc8 = buildSensorYearCollection(i, startDay, endDay, 'LC08', aoi, exclude);
    var lc9 = buildSensorYearCollection(i, startDay, endDay, 'LC09', aoi, exclude)
    var tmp = ee.ImageCollection(lt5.merge(le7).merge(lc8).merge(lc9)); 
    if(first === true){
      var all = tmp;
      first = false;
    } else{
      all = all.merge(tmp);
    }
  }
  
  return ee.Dictionary({
    'idList':all.toList(all.size().add(1)).map(getImgID),
    'collection':all
  });
};
exports.getCollectionIDlist = getCollectionIDlist;



//------ FUNCTION TO COUNT NUMBER OF UNMASKED PIXELS IN AN INTRA ANNUAL COLLECTION ------------------------------
var countClearViewPixels = function(intraAnnualSRcollection){
  var binary = intraAnnualSRcollection.map(function(img){
    return img.select(0)
              .multiply(0)
              .add(1)
              .unmask(0);
  });
  return binary.sum();
};
exports.countClearViewPixels = countClearViewPixels;

//------ FUNCTION TO BUILD ANNUAL COLLECTION OF NUMBER OF UNMASKED PIXELS AVAILABLE TO BUILD COMPOSITE ------------------------------
var buildClearPixelCountCollection = function(startYear, endYear, startDay, endDay, aoi, maskThese) {
  var dummyCollection = ee.ImageCollection([ee.Image([0,0,0,0,0,0]).mask(ee.Image(0))]);
  var imgs = [];     
  for (var i = startYear; i <= endYear; i++) {
    var collection = getCombinedSRcollection(i, startDay, endDay, aoi, maskThese, maskThese);
    var imageCount = collection.toList(1).length();
    var finalCollection = ee.ImageCollection(ee.Algorithms.If(imageCount.gt(0), collection, dummyCollection)); 
    var notMaskCount = countClearViewPixels(finalCollection);
    imgs = imgs.concat(notMaskCount.set('system:time_start', (new Date(i,8,1)).valueOf()));      
  }
  return ee.ImageCollection(imgs);
};
exports.buildClearPixelCountCollection = buildClearPixelCountCollection;



var removeImages = function(collection, exclude){
  // ["LANDSAT/LC08/C01/T1_SR/LC08_046028_20170815"](system:id) or [LC08_046028_20170815](system:index)
  // could not get (system:id) to work though, so getting via string split and slice
  if('exclude' in exclude){
    //print('in exclude')
    exclude = exclude.exclude;
    if('imgIds' in exclude){
      //print('in imgIds')
      var excludeList = exclude.imgIds;
      for(var i=0; i<excludeList.length; i++){
        //print('img blah blah')
        collection = collection.filter(ee.Filter.neq('system:index', excludeList[i].split('/').slice(-1).toString())); //system:id
      }
    }
    if('slcOff' in exclude){
      //print('in slcOff')
      if(exclude.slcOff === true){
        //print('slcOff is true')        
        //'SATELLITE' changed to SPACECRAFT_ID and 'SENSING_TIME' to 'SCENE_CENTER_TIME' in collection 2
        collection = collection.filter(ee.Filter.and(ee.Filter.eq('SPACECRAFT_ID', 'LANDSAT_7'), ee.Filter.gt('SCENE_CENTER_TIME', '2003-06-01T00:00')).not());
      }
    }
  }
  return collection;
};
exports.removeImages = removeImages;




// #######################################################################################
// ###### INDEX CALCULATION FUNCTIONS ####################################################
// #######################################################################################

//LAI collecttion harmonized MODIS and AVHRR annual composite at 5.56 km scale

function qualityMaskModis(image) {
  var qa = image.select('FparLai_QC')

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var  MODLAND = 1 << 0;
  var CLOUD = 1 << 3;
  var SCF = 3 << 5

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(MODLAND).eq(0).and(
            qa.bitwiseAnd(CLOUD).eq(0)).and(
            qa.bitwiseAnd(SCF).eq(0))

  // Return the masked and scaled data, without the QA bands.
  return image.updateMask(mask)
  .select(['Lai', 'Fpar']).rename(['LAI', 'FAPAR'])
      .copyProperties(image, ["system:time_start"])
}

function qualityMaskAVHRR(image) {
  var qa = image.select('QA')

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var  cloud = 1 << 1;
  var invalid = 2 << 1;
  var oor = 3 << 1

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloud).eq(0).and(
            qa.bitwiseAnd(invalid).eq(0)).and(
            qa.bitwiseAnd(oor).eq(0))

  // Return the masked and scaled data, without the QA bands.
  return image.updateMask(mask)
  .select(['LAI','FAPAR'])
      .copyProperties(image, ["system:time_start"])
}




function LAIcol(startYear, endYear, startDay, endDay, aoi){
  var proj = ee.ImageCollection("NOAA/CDR/AVHRR/LAI_FAPAR/V4").first().projection()
  var scale = proj.nominalScale()
  /////Paste all functions
  
  var downres = function(img){
    var reproj = img.setDefaultProjection(proj)
    .reduceResolution({
      reducer:ee.Reducer.median() 
    }).reproject({
      crs:proj,
      scale:scale
    })
    return reproj
  }
  
  var filterCollection = function(year, startDay, endDay, aoi){
  return ee.ImageCollection("NOAA/CDR/AVHRR/LAI_FAPAR/V4")
           .filterBounds(aoi)
           .filterDate(year+'-'+startDay, year+'-'+endDay)
           .map(qualityMaskAVHRR)
           
};
  
  
  var buildSensorYearCollection = function(year, startDay, endDay,aoi){
    var startMonth = parseInt(startDay.substring(0, 2));
    var endMonth = parseInt(endDay.substring(0, 2));
    var srCollection;
    if(startMonth > endMonth){
      var oldYear = (parseInt(year)-1).toString();
      var newYear = year;
      var oldYearStartDay = startDay;
      var oldYearEndDay = '12-31';
      var newYearStartDay = '01-01';
      var newYearEndDay = endDay;
      
      var oldYearCollection = filterCollection(oldYear, oldYearStartDay, oldYearEndDay, aoi);
      var newYearCollection = filterCollection(newYear, newYearStartDay, newYearEndDay, aoi);
      
      srCollection = ee.ImageCollection(oldYearCollection.merge(newYearCollection));
    } else {
      srCollection = filterCollection(year, startDay, endDay, aoi);
    }
  
    
    
    return srCollection;
  };
  
  var years_old = ee.List.sequence(startYear, 2005).getInfo()
  
  var LAI = ee.ImageCollection(years_old.map(function(year){
    var image = buildSensorYearCollection(year, startDay, endDay, aoi).median()
    return image.clip(aoi).set('system:time_start', new Date(year, 8, 1).getTime())}))
    .map(downres)
  
  
  var filterCollection = function(year, startDay, endDay, aoi){
    return ee.ImageCollection("MODIS/061/MCD15A3H")
             .filterBounds(aoi)
             .filterDate(year+'-'+startDay, year+'-'+endDay)
             .map(qualityMaskModis)
             
  };
  
  var years_new = ee.List.sequence(2005, endYear).getInfo()
  
  
  var LAI_new = ee.ImageCollection(years_new.map(function(year){
    var image = buildSensorYearCollection(year, startDay, endDay, aoi).median()
    return image.clip(aoi)
    .set('system:time_start', new Date(year, 8, 1).getTime())
    })).map(downres)
  
  
  var LAI = LAI.map(function(image){
             var LAI = image.select('LAI').multiply(0.0001).multiply(0.3025988021153772).add(0.14168120305862122)
             var FAPAR = image.select('FAPAR').multiply(0.001).multiply(0.44193165531725004).add(0.2805638310366593)
             var scaled = ee.Image().addBands([LAI,FAPAR], ['LAI','FAPAR'])
             return scaled.copyProperties(image, ['system:time_start'])
           })
           
  var LAI_new = LAI_new.map(function(image){
             var LAI = image.select('LAI').multiply(0.01).multiply(0.3025988021153772).add(0.14168120305862122)
             var FAPAR = image.select('FAPAR').multiply(0.01).multiply(0.44193165531725004).add(0.2805638310366593)
             var scaled = ee.Image().addBands([LAI,FAPAR], ['LAI','FAPAR'])
             return scaled.copyProperties(image, ['system:time_start'])
           })
  
  
  var merged = LAI.merge(LAI_new)
  return merged.select(['LAI','FAPAR'])
}

exports.LAIcol = LAIcol;


// TASSELLED CAP
var tcTransform = function(img){ 
  var b = ee.Image(img).select(["B1", "B2", "B3", "B4", "B5", "B7"]); // select the image bands
  var brt_coeffs = ee.Image.constant([0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303]); // set brt coeffs - make an image object from a list of values - each of list element represents a band
  var grn_coeffs = ee.Image.constant([-0.1603, -0.2819, -0.4934, 0.7940, -0.0002, -0.1446]); // set grn coeffs - make an image object from a list of values - each of list element represents a band
  var wet_coeffs = ee.Image.constant([0.0315, 0.2021, 0.3102, 0.1594, -0.6806, -0.6109]); // set wet coeffs - make an image object from a list of values - each of list element represents a band
  
  var sum = ee.Reducer.sum(); // create a sum reducer to be applyed in the next steps of summing the TC-coef-weighted bands
  var brightness = b.multiply(brt_coeffs).reduce(sum); // multiply the image bands by the brt coef and then sum the bands
  var greenness = b.multiply(grn_coeffs).reduce(sum); // multiply the image bands by the grn coef and then sum the bands
  var wetness = b.multiply(wet_coeffs).reduce(sum); // multiply the image bands by the wet coef and then sum the bands
  var angle = (greenness.divide(brightness)).atan().multiply(180/Math.PI).multiply(100);
  var tc = brightness.addBands(greenness)
                     .addBands(wetness)
                     .addBands(angle)
                     .select([0,1,2,3], ['TCB','TCG','TCW','TCA']) //stack TCG and TCW behind TCB with .addBands, use select() to name the bands
                     .set('system:time_start', img.get('system:time_start'));
  return tc;
};

// NBR
var nbrTransform = function(img) {
    var nbr = img.normalizedDifference(['B4', 'B7']) // calculate normalized difference of B4 and B7. orig was flipped: ['B7', 'B4']
                 .multiply(1000) // scale results by 1000
                 .select([0], ['NBR']) // name the band
                 .set('system:time_start', img.get('system:time_start'));
    return nbr;
};

// NDFI - from CODED utility (original: users/bullocke/coded:coded/miscUtilities)
var ndfiTransform = function(img) {

  // pre-defined endmembers
  var params = ee.Dictionary({
    'cfThreshold': 0.01, // CLOUD THRESHOLD 
    'soil': [2000, 3000, 3400, 5800, 6000, 5800],
    'gv': [500, 900, 400, 6100, 3000, 1000],
    'npv': [1400, 1700, 2200, 3000, 5500, 3000],
    'shade': [0, 0, 0, 0, 0, 0],
    'cloud': [9000, 9600, 8000, 7800, 7200, 6500]
    });
    
    /* Utility function for calculating spectral indices */
    var gv = params.get('gv');
    var shade = params.get('shade');
    var npv = params.get('npv');
    var soil = params.get('soil');
    var cloud = params.get('cloud');
    //var cfThreshold = ee.Image.constant(params.get('cfThreshold'))
    /*  Do spectral unmixing on a single image  */
    var unmixImage = ee.Image(img).unmix([gv, shade, npv, soil, cloud], true,true)
                    .rename(['band_0', 'band_1', 'band_2','band_3','band_4']);
    var newImage = ee.Image(img).addBands(unmixImage);
    //var mask = newImage.select('band_4').lt(cfThreshold)
  
    var ndfi = unmixImage.expression(
      '((GV / (1 - SHADE)) - (NPV + SOIL)) / ((GV / (1 - SHADE)) + NPV + SOIL)', {
        'GV': unmixImage.select('band_0'),
        'SHADE': unmixImage.select('band_1'),
        'NPV': unmixImage.select('band_2'),
        'SOIL': unmixImage.select('band_3')
      }); 
    var ndvi = ee.Image(img).normalizedDifference(['B4','B3']).rename('NDVI')
    var evi = ee.Image(img).expression(
          'float(2.5*(((B4/10000) - (B3/10000)) / ((B4/10000) + (6 * (B3/10000)) - (7.5 * (B1/10000)) + 1)))',
          {
              'B4': ee.Image(img).select(['B4']),
              'B3': ee.Image(img).select(['B3']),
              'B1': ee.Image(img).select(['B1'])
          }).rename('EVI');    
          
    var toExp = newImage
          .addBands([ndfi.rename(['NDFI']), ndvi, evi])
          .select(['band_0','band_1','band_2','band_3','NDFI','NDVI','EVI','B1','B2','B3','B4','B5'])
          .rename(['GV','Shade','NPV','Soil','NDFI','NDVI','EVI','Blue','Green','Red','NIR','SWIR1']); 
          //.updateMask(mask)

  toExp = toExp.select(['NDFI'])
              .multiply(1000)
              .set('system:time_start', img.get('system:time_start'));
  return toExp;

};

// NDVI
var ndviTransform = function(img){ 
  var ndvi = img.normalizedDifference(['B4', 'B3']) // calculate normalized dif between band 4 and band 3 (B4-B3/B4_B3)
                .multiply(1000) // scale results by 1000
                .select([0], ['NDVI']) // name the band
                .set('system:time_start', img.get('system:time_start'));
  return ndvi;
};


// NDSI
var ndsiTransform = function(img){ 
  var ndsi = img.normalizedDifference(['B2', 'B5']) // calculate normalized dif between band 2 and band 5 (B2-B5/B2+B5)
                .multiply(1000) // scale results by 1000
                .select([0], ['NDSI']) // name the band
                .set('system:time_start', img.get('system:time_start'));
  return ndsi;
};

// NDMI (sometimes referred to as NDII, normalized difference infrared index)
var ndmiTransform = function(img) {
    var ndmi = img.normalizedDifference(['B4', 'B5']) // calculate normalized difference of B4 and B7. orig was flipped: ['B7', 'B4']
                 .multiply(1000) // scale results by 1000
                 .select([0], ['NDMI']) // name the band
                 .set('system:time_start', img.get('system:time_start'));
    return ndmi;
};

// EVI
var eviTransform = function(img) {
  var evi = img.expression(
      '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
        'NIR': img.select('B4'),
        'RED': img.select('B3'),
        'BLUE': img.select('B1')
  })
  .multiply(1000) // scale results by 1000
  .select([0], ['EVI']) // name the band
  .set('system:time_start', img.get('system:time_start')); 
  return evi;
};

// GNDVI
var gndviTransform = function(img){ 
  var gndvi = img.normalizedDifference(['B4', 'B2']) // calculate normalized dif between band 4 and band 3 (B4-B3/B4_B3)
                .multiply(1000) // scale results by 1000
                .select([0], ['GNDVI']) // name the band
                .set('system:time_start', img.get('system:time_start'));
  return gndvi;
};

// CALCULATE A GIVEN INDEX
var calcIndex = function(img, index, flip){
  // make sure index string in upper case
  index = index.toUpperCase();
  
  // figure out if we need to calc tc
  var tcList = ['TCB', 'TCG', 'TCW', 'TCA'];
  var doTC = tcList.indexOf(index);
  if(doTC >= 0){
    var tc = tcTransform(img);
  }
  
  // need to flip some indices if this is intended for segmentation
  var indexFlip = 1;
  if(flip == 1){
    indexFlip = -1;
  }
  
  // need to cast raw bands to float to make sure that we don't get errors regarding incompatible bands
  // ...derivations are already float because of division or multiplying by decimal
  var indexImg;
  switch (index){
    case 'B1':
      indexImg = img.select(['B1']).float();//.multiply(indexFlip);
      break;
    case 'B2':
      indexImg = img.select(['B2']).float();//.multiply(indexFlip);
      break;
    case 'B3':
      indexImg = img.select(['B3']).float();//.multiply(indexFlip);
      break;
    case 'B4':
      indexImg = img.select(['B4']).multiply(indexFlip).float();
      break;
    case 'B5':
      indexImg = img.select(['B5']).float();//.multiply(indexFlip);
      break;
    case 'B7':
      indexImg = img.select(['B7']).float();//.multiply(indexFlip);
      break;
    case 'NBR':
      indexImg = nbrTransform(img).multiply(indexFlip);
      break;
    case 'NDMI':
      indexImg = ndmiTransform(img).multiply(indexFlip);
      break;
    case 'NDVI':
      indexImg = ndviTransform(img).multiply(indexFlip);
      break;
    case 'NDSI':
      indexImg = ndsiTransform(img).multiply(indexFlip);
      break;
    case 'EVI':
      indexImg = eviTransform(img).multiply(indexFlip);
      break;
    case 'GNDVI':
      indexImg = gndviTransform(img).multiply(indexFlip);
      break;
    case 'TCB':
      indexImg = tc.select(['TCB'])//.multiply(indexFlip);
      break;
    case 'TCG':
      indexImg = tc.select(['TCG']).multiply(indexFlip);
      break;
    case 'TCW':
      indexImg = tc.select(['TCW']).multiply(indexFlip);
      break;
    case 'TCA':
      indexImg = tc.select(['TCA']).multiply(indexFlip);
      break;
     case 'NDFI':
      indexImg = ndfiTransform(img).multiply(indexFlip);
      break;
    default:
      print('The index you provided is not supported');
  };
  

  return indexImg.set('system:time_start', img.get('system:time_start'));
};


exports.calcIndex = calcIndex;



var standardize = function(collection){
  var mean = collection.reduce(ee.Reducer.mean());
  var stdDev = collection.reduce(ee.Reducer.stdDev());
  
  var meanAdj = collection.map(function(img){
    return img.subtract(mean).set('system:time_start', img.get('system:time_start'));
  });
  
  return meanAdj.map(function(img){
    return img.divide(stdDev).set('system:time_start', img.get('system:time_start'));
  });
};

exports.standardize = standardize;

var standardize = function(collection){
  var mean = collection.reduce(ee.Reducer.mean());
  var stdDev = collection.reduce(ee.Reducer.stdDev());
  
  var meanAdj = collection.map(function(img){
    return img.subtract(mean).set('system:time_start', img.get('system:time_start'));
  });
  
  return meanAdj.map(function(img){
    return img.divide(stdDev).set('system:time_start', img.get('system:time_start'));
  });
};

exports.standardize = standardize;



// TRANSFORM AN ANNUAL SR COLLECTION TO AN ANNUAL COLLECTION OF SELECTED INDICES OR BANDS
var transformSRcollection = function(srCollection, bandList){
  return srCollection.map(function(img){
    var allStack = calcIndex(img, bandList[0], 0);
    for(var band=1; band < bandList.length; band++){
      var bandImg = calcIndex(img, bandList[band], 0);
      allStack = allStack.addBands(bandImg);
    }
    return allStack.set('system:time_start', img.get('system:time_start'));
  });
};

exports.transformSRcollection = transformSRcollection;



//##################################################################################################################
//##
//##################################################################################################################


function createTrainingImage(yearPassed, dataset, trainingClassLevel, aoi, customClassLevels){
  // This function takes a land cover map, remaps the classes to a required level and outputs a training image
  // year: String, Classes: int: 0-3
  // Custom class levels allows the user to define a class remapping scheme that acheives separation in classes of interest
  // TODO: condition based on yearPassed which image is (year) used 
  
  // if(yearPassed <= 1990){
  //   var yearUsed=yearPassed
  // }
  
  var yearUsed=yearPassed;
    
  if (dataset == 'CORINE'){
    var originalClasses = [111, 112, 121, 122, 123, 124, 131, 132, 133, 141, 142, 211, 212, 213, 221, 222, 223, 231, 241, 242, 243, 244, 311, 312, 313, 321, 322, 323, 324, 331, 332, 333, 334, 335, 411, 412, 421, 422, 423, 511, 512, 521, 522, 523];
    
    
    // Level one classes Corine:
    
    if (trainingClassLevel == 0){
      //The class level 1 takes the original 43 classes down to 5 levels: ['artificial', 'Agricultural', 'Forest and semi-natural', 'waterbodies', 'Wetlands']
      
      var newClasses = [1,1,1,1,1,1,1,1,1,1,1,    //Artificial (11)
                        2,2,2,2,2,2,2,2,2,2,2,    //Agricultural(11)
                        3,3,3,3,3,3,3,3,3,3,3,3,  //Forest and Semi Natural(12)
                        0,0,0,0,0,                //Waterbodies(5)
                        0,0,0,0,0];                //Wetlands(5)
      
      
      
      var trainingImage = ee.Image('COPERNICUS/CORINE/V20/100m/'+yearUsed)
                          .clip(aoi)
                          .remap(originalClasses, newClasses)
                          .rename('landcover');
                          
      var mask = trainingImage.select('landcover').gt(0)
      
      var trainingImage = trainingImage.updateMask(mask)
      
    }
    
    
    // Level one classes Corine:
    
    if (trainingClassLevel == 1){
      //The class level 1 takes the original 43 classes down to 5 levels: ['artificial', 'Agricultural', 'Forest and semi-natural', 'waterbodies', 'Wetlands']
      
      var newClasses = [1,1,1,1,1,1,1,1,1,1,1,    //Artificial (11)
                        2,2,2,2,2,2,2,2,2,2,2,    //Agricultural(11)
                        3,3,3,3,3,3,3,3,3,3,3,3,  //Forest and Semi Natural(12)
                        4,4,4,4,4,                //Waterbodies(5)
                        5,5,5,5,5];                //Wetlands(5)
                              
      
      var trainingImage = ee.Image('COPERNICUS/CORINE/V20/100m/'+yearUsed)
                          .clip(aoi)
                          .remap(originalClasses, newClasses)
                          .rename('landcover');
    }
    
    // Level two classes Corine:
    
    if (trainingClassLevel == 2){
      //The class level 2 takes the original 44 classes down to 15 levels as described below
      
      var newClasses = [1,1,          //Urban Fabric (2)
                        2,2,2,2,      //Industrial, Commercial & TransportUnits (4)
                        3,3,3,        //Mine, dump & construction (3)
                        4,4,          //Artificial, NonAgricultural, vegetated (2)
                        5,5,5,        //Arable land (3)
                        6,6,6,        //Permanent crops (3)
                        7,            //Pastures (1)
                        8,8,8,8,      //Heterogeneous Agricultural Areas (4)
                        9,9,9,        //Forests (3)
                        10,10,10,10,  //Scrub and/or herbaceous vegetation associations (4)
                        11,11,11,11,11,//Open Spaces with litte or no vegetation (5)
                        12,12,        //Inland wetlands (2)
                        13,13,13,     //Maritime wetlands (3)
                        14,14,        //Inland waters (2)
                        15,15,15] ;    //Marine waters (3)
      
      var trainingImage = ee.Image('COPERNICUS/CORINE/V20/100m/'+yearUsed)
                          .clip(aoi)
                          .remap(originalClasses, newClasses)
                          .rename('landcover');
    }
    
    //Level three classes, original 44
    
    if (trainingClassLevel == 3){
      //The class level 3 returns the original 44 classes
      
      var trainingImage = ee.Image('COPERNICUS/CORINE/V20/100m/'+yearUsed)
                          .clip(aoi);
    }
    
    //Level four (customClassLevels is a list passed by the user, of length 44)
    if (trainingClassLevel == 4){
      var newClasses = customClassLevels;
      var trainingImage = ee.Image('COPERNICUS/CORINE/V20/100m/'+yearUsed)
                          .clip(aoi)
                          .remap(originalClasses, newClasses)
                          .rename('landcover');
    }
  }
  
  //The clipped and remapped image is returned. Working well with all 4 conditions. 
  return trainingImage
}



exports.createTrainingImage = createTrainingImage;

function addTerrainBand(image, aoi){

  var DEM = ee.Image("USGS/SRTMGL1_003").clip(aoi);
  var elev = DEM.select('elevation').rename('elev');
  var slop = ee.Terrain.slope(DEM.select('elevation')).rename('slope');
  var add = image.addBands(elev).addBands(slop);
  
  return add
  
};


exports.addTerrainBand = addTerrainBand;

function genGCP(trainingClassImage, imageToClassify, numClasses, split, tileScale, aoi, dist, pct){

  if (numClasses == 3){
    var sample = trainingClassImage.stratifiedSample({
                    numPoints: ee.Number(5000).round(), 
                    classBand: "landcover",
                    region:aoi, 
                    scale: 30,
                    projection:trainingClassImage.projection(),
                    classValues:[1, 2, 3],
                    classPoints:[ee.Number(5000*0.15).round(), ee.Number(5000*0.25).round(),ee.Number(5000*0.60).round()],
                    geometries:true
                  });
  }
  // Extract random samples based on the 'landcover' band  
  // Most commnoly a memory limit is reached when the number of points exceeds 5000
  if (numClasses == 5 && dist == 'weighted'){
    var sample = trainingClassImage.stratifiedSample({
                    numPoints: ee.Number(5000).round(), 
                    classBand: "landcover",
                    region:aoi, 
                    scale: 30,
                    projection:trainingClassImage.projection(),
                    classValues:[1, 2, 3, 4, 5],
                    classPoints:[ee.Number(5000*pct[0]).round(), ee.Number(5000*pct[1]).round(),ee.Number(5000*pct[2]).round(),ee.Number(5000*pct[3]).round(), ee.Number(5000*pct[4]).round()],
                    geometries:true
                  });
  }
  
  if (numClasses == 5 && dist == 'balanced'){
    var sample = trainingClassImage.stratifiedSample({
                    numPoints: ee.Number((5000/numClasses)).round(), 
                    classBand: "landcover",
                    region:aoi, 
                    scale: 30,
                    projection:trainingClassImage.projection(),
                    geometries:true
                  });
  }
  
  
  var addRand = sample.randomColumn('random');
  
  var trainingGcp = addRand.filter(ee.Filter.lt('random', split));
  var testingGcp = addRand.filter(ee.Filter.gte('random', split));
  
  
  var training = imageToClassify.sampleRegions({
    collection: trainingGcp,
    properties: ['landcover'],
    scale: 10,
    tileScale:1,
  });
  
  
  var test= imageToClassify.sampleRegions({
    collection: testingGcp,
    properties: ['landcover'],
    scale:10,
    tileScale:1,
  });

  return {training:training, testing:test}
  
}

exports.genGCP = genGCP;

function classifier(toClassify, trainingGcp, tuned, params){
  // Classify function currently usese the Copernicus Corine Dataset. This is writtent to investigate Decae-over-decade change in line with Discharge data available
  // The classifier takes a year that corresponds to the year 
  // Optimal nTrees (approx = 70) determined in separate testing
    
    if (tuned == 'not_tuned'){
    var classifier = ee.Classifier.smileRandomForest(50)
                                  .train({
                                    features:trainingGcp, 
                                    classProperty:'landcover', 
                                    inputProperties: toClassify.bandNames(),
                                    });
    }
    
    if (tuned == 'tuned'){
        var classifier = ee.Classifier.smileRandomForest({numberOfTrees:params[0],
                                                     variablesPerSplit:params[1],
                                                     minLeafPopulation:params[2], 
                                                     bagFraction:params[3],
                                                     maxNodes:params[4],
                                                     seed:params[5]
                                                         })
                                  .train({
                                    features:trainingGcp, 
                                    classProperty:'landcover', 
                                    inputProperties: toClassify.bandNames(),
                                    });
    }
  
  
  
  return classifier
  
};

exports.classifier = classifier;


function classArea(classImage, scale, geom){ 
  var areaImage = ee.Image.pixelArea().addBands(classImage)
  var areas = areaImage.reduceRegion({
        reducer: ee.Reducer.sum().group({
        groupField: 1,
        groupName: 'class',
      }),
      geometry: geom,
      scale: scale,
      maxPixels: 1e10
      }); 
  
  var classAreas = ee.List(areas.get('groups'))
  
  var classAreaLists = classAreas.map(function(item) {
    var areaDict = ee.Dictionary(item)
    var classNumber = ee.Number(areaDict.get('class')).format()
    var area = ee.Number(areaDict.get('sum')).divide(1e6)//.round()
    return ee.List([classNumber, area])
  })
  
  var result = ee.Dictionary(classAreaLists.flatten())
  return result
};

exports.classArea = classArea;




function imcolFromAsset(startYear, endYear, aoi, dataset, asset_folder, file_suffix, SMN){
  var imageCollection = ee.ImageCollection([]);
  for (var year = startYear; year <= endYear; year++) {
    var asset_id = '';
    if (dataset == 'Meuse') {
      asset_id = asset_folder + '/'+dataset+'_' + year.toString() + file_suffix;
    } else if (dataset == 'CAMELS_GB') {
      asset_id = asset_folder + '/GB_' + SMN + '_HC_tuned_' + year.toString();
    }
    var image = ee.Image(asset_id).clip(aoi).set('system:time_start', ee.Date.fromYMD(year, 8, 31));
    var image_geom = image.set('system:footprint', aoi);
    imageCollection = imageCollection.merge(image_geom);
  }
  return imageCollection;
}

exports.imcolFromAsset = imcolFromAsset;


function imcolFromAssetHILDA(startYear, endYear, aoi, dataset, asset_folder) {
  var imageCollection = ee.ImageCollection([]);
  
  for (var year = startYear; year <= endYear; year++) {
    var asset_id = asset_folder + '/' + 'clipped_HILDAplus_geotiff_' + dataset + '_' + year.toString() + '_4326';
    print(asset_id)
    var date = ee.Date.fromYMD(year, 8, 31).advance(22, 'hour');
    var timeStart = date.millis();
    var image = ee.Image(asset_id).clip(aoi).set('system:time_start', timeStart);
    var image_geom = image.set('system:footprint', aoi);
    
    imageCollection = imageCollection.merge(image_geom);
  }
  
  return imageCollection;
}

exports.imcolFromAssetHILDA = imcolFromAsset;





var unionImages = function(image) {
  // Get the primary and secondary images
  var primary = ee.Image(image.get('primary'));
  var secondary = ee.Image(image.get('secondary'));

  // Get the union of the footprints
  var unionedFootprint = primary.geometry().union(secondary.geometry());

  // Clip the images to the unioned footprint
  var primaryClipped = primary.clip(unionedFootprint);
  var secondaryClipped = secondary.clip(unionedFootprint);

  // Unmask the images
  primaryClipped = primaryClipped.unmask(0);
  secondaryClipped = secondaryClipped.unmask(0);

  // Merge the images
  var merged = primaryClipped.add(secondaryClipped);

  // Set the system:footprint property of the merged image
  merged = merged.set('system:footprint', unionedFootprint).set('system:time_start', ee.Date(primary.get('system:time_start')).millis());

  return merged;
};

function unionCollections(imageCollection_a, imageCollection_b){
    var filter = ee.Filter.equals({
                                  leftField: 'system:time_start',
                                  rightField: 'system:time_start'
                                });

    // Create the join
    var simpleJoin = ee.Join.inner();
    var innerJoin = ee.ImageCollection(simpleJoin.apply(imageCollection_a, imageCollection_b, filter));

    // Map over the joined collection to union the images
    var unioned = innerJoin.map(unionImages);
    
    return unioned
}

exports.unionCollections = unionCollections;