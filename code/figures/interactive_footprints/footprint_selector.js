var display = data_display.data;
var ex_display = exshift_display.data;

var source_data = data.data;
var exshift_data = exshift.data;

var prom = prom_selector.value;
var gc = gc_selector.value;
var replicate = rep_selector.value;
var d = d_selector.value;


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



// Custom Function Definitions
function getAllIndexes(arr, val) {
    var indices = [], i = -1;
    while ((i = arr.indexOf(val, i+1)) != -1){
        indices.push(i);
    }
    return indices;
}