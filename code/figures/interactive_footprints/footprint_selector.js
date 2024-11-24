var display = data_display.data;
var ex_display = exshift_display.data;
var rep_comparer = replicate_comparer.data;
var rep_comparer_line = replicate_comparer_line.data;
var ang_display = angle_display.data;


var source_data = data.data;
var exshift_data = exshift.data;
var shuffle_data = shuffle.data;

var prom = prom_selector.value;
var gc = gc_selector.value;
var replicate = rep_selector.value;
var d = Number(d_selector.value);
var background_sub = background_selector.value;
var sigma_multiplier = sigma_selector.value;


// Get data for footprints
var prom_inds = getAllIndexes(source_data['promoter'], prom);
var gc_inds = getAllIndexes(source_data['growth_condition'], gc);
var rep_inds = getAllIndexes(source_data['replicate'], replicate);

var filteredArray = prom_inds.filter(value => gc_inds.includes(value));
var filteredArray = filteredArray.filter(value => rep_inds.includes(value));

display['pos'] = source_data['pos'][filteredArray];
display['mut_info'] = source_data['mut_info'][filteredArray];

// Get shuffled data
var prom_inds = getAllIndexes(shuffle_data['promoter'], prom);
var gc_inds = getAllIndexes(shuffle_data['gc_name'], gc);
var rep_inds = getAllIndexes(shuffle_data['rep'], replicate);

var filteredArray_shuffle = prom_inds.filter(value => gc_inds.includes(value));
var filteredArray_shuffle = filteredArray.filter(value => rep_inds.includes(value));


if (d == 0) {
    display['pos'] = source_data['pos'][filteredArray];
    display['mut_info'] = source_data['mut_info'][filteredArray];
    if (background_sub == "yes"){
        var mut_info = [];
        for (var i = 0; i < (source_data['mut_info'][filteredArray].length); i=i+1) {
            mut_info.push(Math.max(display['mut_info'][i] - (sigma_multiplier * shuffle_data['sigma'][filteredArray][i] + shuffle_data['mu'][filteredArray][i]), 0))
        }
        display['mut_info'] = mut_info;
    }
  } else {
    var mut_info = [];
    var pos = [];
    for (var i = d; i < (source_data['mut_info'][filteredArray].length - d); i=i+1) {
        pos.push(source_data['pos'][filteredArray][i]);
        mut_info.push(source_data['mut_info'][filteredArray].slice(i-d, i+d+1).reduce((a,b)=>a+b) / d);
    }
    display['pos'] = pos;
    display['mut_info'] = mut_info;
    if (background_sub == "yes"){
        var threshold = [];
        for (var i = 0; i < (shuffle_data['mu'][filteredArray].length); i=i+1) {
            threshold.push(shuffle_data['mu'][filteredArray][i] + shuffle_data['sigma'][filteredArray][i] * sigma_multiplier)
        }
        
        var mut_info_shuffle = [];
        var mut_info_filt = [];
        
        for (var i = d; i < (shuffle_data['mu'][filteredArray].length - d); i=i+1) {
            mut_info_shuffle.push(threshold.slice(i-d, i+d+1).reduce((a,b)=>a+b) / d);
        }
        for (var i = 0; i < (mut_info_shuffle.length); i=i+1) {
            mut_info_filt.push(Math.max(display['mut_info'][i] - mut_info_shuffle[i], 0))
        }
        console.log(mut_info_shuffle)
        display['mut_info'] = mut_info_filt;
    }
}


// data for replicate plot
var prom_inds = getAllIndexes(source_data['promoter'], prom);
var gc_inds = getAllIndexes(source_data['growth_condition'], gc);
//var d_inds = getAllIndexes(source_data['d'], d);

var rep_inds_1 = getAllIndexes(source_data['replicate'], '1');
var rep_inds_2 = getAllIndexes(source_data['replicate'], '2');

var filteredArray = prom_inds.filter(value => gc_inds.includes(value));
//var filteredArray = filteredArray.filter(value => d_inds.includes(value));

var filteredArray_1 = filteredArray.filter(value => rep_inds_1.includes(value));
var filteredArray_2 = filteredArray.filter(value => rep_inds_2.includes(value));

if (d == 0) {
    var rep1 = source_data['mut_info'][filteredArray_1];
    var pos1 = source_data['pos'][filteredArray_1];
    var rep2 = source_data['mut_info'][filteredArray_2];
    var pos2 = source_data['pos'][filteredArray_2];
  } else {
    var rep1 = [];
    var pos1 = [];
    var rep2 = [];
    var pos2 = [];
    for (var i = d; i < (source_data['pos'][filteredArray_1].length - d); i=i+1) {
        pos1.push(source_data['pos'][filteredArray_1][i]);
        pos2.push(source_data['pos'][filteredArray_2][i]);
        rep1.push(source_data['mut_info'][filteredArray_1].slice(i-d, i+d+1).reduce((a,b)=>a+b) / d);
        rep2.push(source_data['mut_info'][filteredArray_2].slice(i-d, i+d+1).reduce((a,b)=>a+b) / d);
    }
}

var rep1_sum = rep1.reduce((partialSum, a) => partialSum + a, 0);
var rep2_sum = rep2.reduce((partialSum, a) => partialSum + a, 0);

var rep1_norm = rep1.map(function(item) { return item/rep1_sum } );
var rep2_norm = rep2.map(function(item) { return item/rep2_sum } );

rep_comparer['rep1'] = rep1_norm;
rep_comparer['rep2'] = rep2_norm;


rep_comparer['x'] = addvector(rep1_norm, rep2_norm).map(function(item) { return item/2 } ) ;
rep_comparer['y'] = divvector(subvector(rep1_norm, rep2_norm), addvector(rep1_norm, rep2_norm)).map(function(item) { return item * 2 } );

rep_comparer['pos'] = pos1;

rep_comparer_line['x'] = [0, getMaxOfArray([getMaxOfArray(rep1), getMaxOfArray(rep2)])];
rep_comparer_line['y'] = [0, getMaxOfArray([getMaxOfArray(rep1), getMaxOfArray(rep2)])];


// Get data for expression shifts
var prom_inds = getAllIndexes(exshift_data['promoter'], prom)
var gc_inds = getAllIndexes(exshift_data['growth_condition'], gc)
var rep_inds = getAllIndexes(exshift_data['replicate'], replicate)

var filteredArray = prom_inds.filter(value => gc_inds.includes(value))
var filteredArray = filteredArray.filter(value => rep_inds.includes(value))

var pos = exshift_data['pos'][filteredArray];
var base = exshift_data['base'][filteredArray];
var wt_base = exshift_data['wt_base'][filteredArray];
var expression_shift = exshift_data['expression_shift'][filteredArray];
ex_display['pos'] = pos;
ex_display['base'] = base;
ex_display['wt_base'] = wt_base;
ex_display['expression_shift'] = expression_shift;


//
var prom_inds = getAllIndexes(exshift_data['promoter'], prom)
var gc_inds = getAllIndexes(exshift_data['growth_condition'], gc)
var rep_inds_1 = getAllIndexes(exshift_data['replicate'], '1')
var rep_inds_2 = getAllIndexes(exshift_data['replicate'], '2')

var filteredArray = prom_inds.filter(value => gc_inds.includes(value))
var filteredArray_1 = filteredArray.filter(value => rep_inds_1.includes(value))
var filteredArray_2 = filteredArray.filter(value => rep_inds_2.includes(value))

var expression_shift_1 = exshift_data['expression_shift'][filteredArray_1];
var expression_shift_2 = exshift_data['expression_shift'][filteredArray_2];

var expression_shift_1_std = getStandardDeviation(expression_shift_1);
var expression_shift_2_std = getStandardDeviation(expression_shift_2);

var expression_shift_1_norm = expression_shift_1.map(function(item) { return item/expression_shift_1_std } );
var expression_shift_2_norm = expression_shift_2.map(function(item) { return item/expression_shift_2_std } );

var cos_list = [];
var norms = [];
var ecdf_y = [];

for (var i = 0; i < expression_shift_2_norm.length; i=i+4) {
    var norm_1 = 0;
    var norm_2 = 0;
    var dot_prod = 0;
    for (var j = 0; j < 4; j=j+1){
        norm_1 += expression_shift_1_norm[i+j]**2;
        norm_2 += expression_shift_2_norm[i+j]**2;
        dot_prod += expression_shift_1_norm[i+j] * expression_shift_2_norm[i+j]
    }
    norms.push(Math.sqrt(norm_1 * norm_2));
    cos_list.push(Math.abs(dot_prod) / Math.sqrt(norm_1 * norm_2));
    ecdf_y.push(i/ expression_shift_2_norm.length);
}

ang_display['x'] = norms;
ang_display['y'] = cos_list;
ang_display['ecdf_x'] = cos_list.sort();
ang_display['ecdf_y'] = ecdf_y;

// get letters for x-axis
const map1 = new Map();

for (var i = 0; i < ex_display['wt_base'].length; i=i+4) {
    var j = (i)/4 - 115;
    map1.set(j, ex_display['wt_base'][i]);
}


x_axis.major_label_overrides = map1;


p.reset.emit();
p.change.emit();


//console.log(display)
data_display.change.emit();

//console.log(ex_display)
exshift_display.change.emit();


//console.log(pos1)
//console.log(pos2)

replicate_comparer.change.emit();
replicate_comparer_line.change.emit();
angle_display.change.emit();


// Custom Function Definitions
function getAllIndexes(arr, val) {
    var indices = [], i = -1;
    while ((i = arr.indexOf(val, i+1)) != -1){
        indices.push(i);
    }
    return indices;
}

function getMaxOfArray(numArray) {
    return Math.max.apply(null, numArray);
  }
  
function addvector(a,b){
    return a.map((e,i) => e + b[i]);
}

function subvector(a,b){
    return a.map((e,i) => Math.abs(e - b[i]));
}

function divvector(a,b){
    return a.map((e,i) => e / b[i]);
}

function getStandardDeviation (array) {
    const n = array.length
    const mean = array.reduce((a, b) => a + b) / n
    return Math.sqrt(array.map(x => Math.pow(x - mean, 2)).reduce((a, b) => a + b) / n)
  }