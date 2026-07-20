var ntl = ee.ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG');
ntl = ntl.filterDate('2014-01-01','2026-04-01');
var extractMonthly = function(img){
  var date = img.date().format('YYYY-MM');
  var ntl_band = img.select('avg_rad');
  var stats = ntl_band.reduceRegions({
    collection: jp,
    reducer: ee.Reducer.mean(),
    scale: 500});
    stats = stats.map(function(f){
    return f.set('month',date);});
  return stats;};
var results = ntl.map(extractMonthly).flatten();
print(results.limit(10));
Export.table.toDrive({collection: results,description: 'Japan_Prefecture_NTL',fileFormat: 'CSV'});