/** Example GEE script for rapid hazard response in Hiroshima, Japan
* Rainfall event: 28 June 2018 – 9 July 2018
* In this example, we define an AOI defined in "Geometry Imports" (see 
* imports.js)
* Code written by Mong-Han Huang, Department of Geology, mhhuang@umd.edu
*/
 
var slope_threshold = .5; // unit: degree
var curv_threshold = -0.005; // unit: m/m^2

//define pre-event stack time period
var PreEventTime_1 = '2015-01-01T23:59'; // format: yyyy-mm-dd-HH:MM
var PreEventTime_2 = '2018-06-29T23:59'; // format: yyyy-mm-dd-HH:MM

//define post-event stack time period
var PostEventTime_1 = '2018-07-09T23:59'; // format: yyyy-mm-dd-HH:MM
var PostEventTime_2 = '2020-05-29T23:59'; // format: yyyy-mm-dd-HH:MM


// filter and mask for Sentinel-2 optical data (not discussed in the manuscript)
function maskS2clouds(image) {
  var qa = image.select('QA60');
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000);
}

// filter and mask for Landsat 8 optical data (not discussed in the manuscript)
function maskL8sr(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  // Get the pixel QA band.
  var qa = image.select('pixel_qa');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
               .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
}

// Map the function and take the median.
// Load Sentinel-2 TOA reflectance data.
var dataset_preEvent = ee.ImageCollection('COPERNICUS/S2')
                  .filterDate(PreEventTime_1,PreEventTime_2)
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10)) //only include images with less than 10% clouds
                  .filterBounds(AOI)
                  .map(maskS2clouds);
                 
var dataset_postEvent = ee.ImageCollection('COPERNICUS/S2')
                  .filterDate(PostEventTime_1, PostEventTime_2)
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10)) //only include images with less than 10% clouds
                  .filterBounds(AOI)
                  .map(maskS2clouds);


// LOAD Shuttle Radar Topography Mission (SRTM) Digital Elevation Model (DEM)
var dataset = ee.Image('USGS/SRTMGL1_003');
var elevation = dataset.select('elevation'); 
var slope = ee.Terrain.slope(elevation); //slope in degrees
var mask_slope = slope.gte(slope_threshold); // slope mask with values 0 or 1
var slope_mask = slope.updateMask(mask_slope);

// Calculate curvature
// Define a Gaussian kernel for smoothing. This step helps reduce noise in the curvature maps
var smooth_curv = ee.Kernel.gaussian({
  radius: 120,
  sigma: 60,
  units: 'meters',
  normalize: true,
});
var xyDemGrad = elevation.convolve(smooth_curv).gradient();
var xGradient = xyDemGrad.select('x').gradient();
var yGradient = xyDemGrad.select('y').gradient();
var curvature = xGradient.select('x').add(yGradient.select('y'));
var mask_curvature = curvature.gte(curv_threshold);
var curvature_mask = slope.updateMask(mask_curvature);

// Define a Gaussian kernel to reduce noise in S1 scenes
var smooth_S1 = ee.Kernel.gaussian({
  radius: 50,
  sigma: 20,
  units: 'meters',
  normalize: true,
});


// LOAD Sentinel-1 (S1) amplitude data VH polarization
var imgVH = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .select('VH')
        .filterBounds(AOI)
        .map(function(image) {
          var edge = image.lt(-30.0); //remove low edge values as suggested by GEE
          var maskedImage = image.mask().and(edge.not());
          return image.updateMask(maskedImage);
        });

var des = imgVH.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING')); //descending acquisition geometry data
var asc = imgVH.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));  //ascending acquisition geometry data

var PreEventPeriod = ee.Filter.date(PreEventTime_1,PreEventTime_2);
var PostEventPeriod = ee.Filter.date(PostEventTime_1,PostEventTime_2);

// median S1 SAR Amplitude
var PreEventPeriod_asc = ee.Image.cat(asc.filter(PreEventPeriod).median());
var PreEventPeriod_des = ee.Image.cat(des.filter(PreEventPeriod).median());

//uncomment to smooth amplitude stack to reduce noise
//var PreEventPeriod_asc = ee.Image.cat(asc.filter(PreEventPeriod).median()).convolve(smooth_S1);
//var PreEventPeriod_des = ee.Image.cat(des.filter(PreEventPeriod).median()).convolve(smooth_S1);

// without amplitude filter
var PostEventPeriod_asc = ee.Image.cat(asc.filter(PostEventPeriod).median());
var PostEventPeriod_des = ee.Image.cat(des.filter(PostEventPeriod).median());

//uncomment to filter amplitude image
//var PostEventPeriod_asc = ee.Image.cat(asc.filter(PostEventPeriod).median()).convolve(smooth_S1);
//var PostEventPeriod_des = ee.Image.cat(des.filter(PostEventPeriod).median()).convolve(smooth_S1);

// print out image information (number of images, image file name)
var num_asc_pre = asc.filter(PreEventPeriod).filterBounds(AOI);
var num_des_pre = des.filter(PreEventPeriod).filterBounds(AOI);
var num_asc_post = asc.filter(PostEventPeriod).filterBounds(AOI);
var num_des_post = des.filter(PostEventPeriod).filterBounds(AOI);

var count_asc_pre = num_asc_pre.sort('system:time_start').toList(5000,0).length(); // 5000 means size of the list
var count_des_pre = num_des_pre.sort('system:time_start').toList(5000,0).length(); // 5000 means size of the list
var count_asc_post = num_asc_post.sort('system:time_start').toList(5000,0).length(); // 5000 means size of the list
var count_des_post = num_des_post.sort('system:time_start').toList(5000,0).length(); // 5000 means size of the list
print("Number of pre-event ascending images: ", count_asc_pre  );
print("Number of pre-event descending images: ", count_des_pre  );
print("Number of post-event ascending images: ", count_asc_post  );
print("Number of post-event ascending images: ", count_des_post  );

print(num_asc_pre,'Pre-event ascending');
print(num_des_pre,'Pre-event descending');
print(num_asc_post,'Post-event ascending');
print(num_des_post,'Post-event descending');


// calculate difference between pre and post event SAR amplitude
var Diff_des = PreEventPeriod_des.subtract(PostEventPeriod_des);
var Diff_asc = PreEventPeriod_asc.subtract(PostEventPeriod_asc);
var Sum_des_asc = (Diff_asc.add(Diff_des)).divide(2);

// define color palette ('ffffff' is white, 'ff0000' is red, '0000ff' is blue, '000000' is black)
var ColorScale = {min: 0, max: 2, palette: ['ffffff','ffffff','ff0000']}; // for amplitude change (A_ratio)
var SlopeColorScale = {min: 0, max: 60, palette: ['ffffff','ff0000','ff0000']}; // for slope
var Grey = {min: -20,max: -9};
var ColorCurv = {min: -0.02, max:0.02, palette: ['ff0000','ffffff','0000ff']}; // for curvature 
var rgbVis = {min: 0.0, max: 0.18, bands: ['B4', 'B3', 'B2']}; // for sentinel-2

///////////////////////////////////////////////////////////////////////////
Map.centerObject(AOI, 14); //zooms to center of AOI after clicking "run". The number determines the zoom level.

//// show Sentinel-2 images
Map.addLayer(dataset_preEvent.median(), rgbVis , 'S2 pre-event',true);
Map.addLayer(dataset_postEvent.median(), rgbVis, 'S2 post-event',true);

//Map.addLayer(PreEventPeriod_asc, Grey, 'Asc - pre amp', true);
//Map.addLayer(PreEventPeriod_des, Grey, 'Des - pre amp', true);
//Map.addLayer((PreEventPeriod_des.add(PreEventPeriod_asc)).divide(2), Grey, 'SUM - pre amp', true);
//Map.addLayer(PostEventPeriod_asc, Grey, 'Asc - post amp', true);
//Map.addLayer(PostEventPeriod_des, Grey, 'Des - post amp', true);
//Map.addLayer((PostEventPeriod_des.add(PostEventPeriod_asc)).divide(2), Grey, 'SUM - post amp', true);
//Map.addLayer(Diff_des, ColorScale, 'SUM: Event Period', true);
//Map.addLayer(Diff_asc, ColorScale, 'SUM: Event Period', true);

//// show DEM masks
//Map.addLayer(curvature, ColorCurv,'SRTM curvature',false);
//Map.addLayer(curvature_mask.updateMask(slope_mask), {min: -1, max:1, palette: ['ff0000','ffffff','0000ff']},'SRTM DEM mask',false);
//Map.addLayer(mask_slope, SlopeColorScale, 'SRTM slope mask', false);

//// show DPM products (with or without a DEM mask)
Map.addLayer(Sum_des_asc, ColorScale, 'S1 SAR amplitude change w/o mask', false);
//Map.addLayer(Sum_des_asc.updateMask(mask_curvature), ColorScale, 'S1 SAR amplitude change curvature mask only', true);
//Map.addLayer(Sum_des_asc.updateMask(mask_slope), ColorScale, 'S1 SAR amplitude change slope mask only', true);
Map.addLayer(Sum_des_asc.updateMask(mask_slope).updateMask(mask_curvature), ColorScale, 'S1 SAR amplitude change w/mask', true);

/////// Export as Geotiff //////////////
// it will save to your google drive as geotiff
Export.image.toDrive({
  image: Sum_des_asc,
  description: 'SAR_amplitude_change',
  scale: 10,
  fileFormat: 'GeoTIFF',
  region: AOI
});


// to comment out line blocks
/*

*/

