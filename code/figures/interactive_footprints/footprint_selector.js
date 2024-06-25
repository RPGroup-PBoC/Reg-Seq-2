var display = data_display.data;
var ex_display = exshift_display.data;
var rep_comparer = replicate_comparer.data;
var rep_comparer_line = replicate_comparer_line.data;

var source_data = data.data;
var exshift_data = exshift.data;

var prom = prom_selector.value;
var gc = gc_selector.value;
var replicate = rep_selector.value;
var d = d_selector.value;


// Get data for footprints
var prom_inds = getAllIndexes(source_data['promoter'], prom);
var gc_inds = getAllIndexes(source_data['growth_condition'], gc);
var rep_inds = getAllIndexes(source_data['replicate'], replicate);
var d_inds = getAllIndexes(source_data['d'], d);

var filteredArray = prom_inds.filter(value => gc_inds.includes(value));
var filteredArray = filteredArray.filter(value => rep_inds.includes(value));
var filteredArray = filteredArray.filter(value => d_inds.includes(value));


var pos = source_data['pos'][filteredArray];
var mut_info = source_data['mut_info'][filteredArray];
display['pos'] = pos;
display['mut_info'] = mut_info;


// data for replicate plot
var prom_inds = getAllIndexes(source_data['promoter'], prom);
var gc_inds = getAllIndexes(source_data['growth_condition'], gc);
var d_inds = getAllIndexes(source_data['d'], d);

var rep_inds_1 = getAllIndexes(source_data['replicate'], '1');
var rep_inds_2 = getAllIndexes(source_data['replicate'], '2');

var filteredArray = prom_inds.filter(value => gc_inds.includes(value));
var filteredArray = filteredArray.filter(value => d_inds.includes(value));

var filteredArray_1 = filteredArray.filter(value => rep_inds_1.includes(value));
var rep1 = source_data['mut_info'][filteredArray_1];
var filteredArray_2 = filteredArray.filter(value => rep_inds_2.includes(value));
var rep2 = source_data['mut_info'][filteredArray_2];

rep_comparer['rep1'] = rep1;
rep_comparer['rep2'] = rep2;
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

//p_exshift.xaxis.major_label_text_font_size = "16pt"
const map1 = new Map();

for (var i = 0; i < ex_display['wt_base'].length; i=i+4) {
    var j = (i)/4 - 115;
    map1.set(j, ex_display['wt_base'][i]);
}


x_axis.major_label_overrides = map1;




p.reset.emit();
p.change.emit();


console.log(display)
data_display.change.emit();

console.log(ex_display)
exshift_display.change.emit();

console.log(rep_comparer)
console.log(rep_comparer_line['x'])

replicate_comparer.change.emit();
replicate_comparer_line.change.emit();



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
  