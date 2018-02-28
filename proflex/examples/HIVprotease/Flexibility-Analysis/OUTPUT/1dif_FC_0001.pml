delete FlexClustObj
delete FC*
delete smallFC*
delete RC*
delete HPHOB3
load 1dif_flex_0001.pdb, FlexClustObj


select FC1, FlexClustObj and id 3070-3086
color wheat, FC1
show sticks, FC1
show cartoon, FC1
cartoon tube, FC1

select FC3, FlexClustObj and id 3093-3110
color lightblue, FC3
show sticks, FC3
show cartoon, FC3
cartoon tube, FC3

select FC4, FlexClustObj and id 3111-3122
color paleyellow, FC4
show sticks, FC4
show cartoon, FC4
cartoon tube, FC4

select FC5, FlexClustObj and id 3123-3185
color lightpink, FC5
show sticks, FC5
show cartoon, FC5
cartoon tube, FC5

select FC7, FlexClustObj and id 3192-3364
color lightorange, FC7
show sticks, FC7
show cartoon, FC7
cartoon tube, FC7

select FC9, FlexClustObj and id 3371-3384
color teal, FC9
show sticks, FC9
show cartoon, FC9
cartoon tube, FC9

select FC10, FlexClustObj and id 3385-3397
color limegreen, FC10
show sticks, FC10
show cartoon, FC10
cartoon tube, FC10

select FC11, FlexClustObj and id 3398-3407
color pink, FC11
show sticks, FC11
show cartoon, FC11
cartoon tube, FC11

select FC13, FlexClustObj and id 3414-3425
color aquamarine, FC13
show sticks, FC13
show cartoon, FC13
cartoon tube, FC13

select FC16, FlexClustObj and id 3438-3480
color tv_orange, FC16
show sticks, FC16
show cartoon, FC16
cartoon tube, FC16

select FC17, FlexClustObj and id 3481-3489
color lime, FC17
show sticks, FC17
show cartoon, FC17
cartoon tube, FC17

select FC18, FlexClustObj and id 3490-3510
color marine, FC18
show sticks, FC18
show cartoon, FC18
cartoon tube, FC18

select FC19, FlexClustObj and id 3511-3519
color splitpea, FC19
show sticks, FC19
show cartoon, FC19
cartoon tube, FC19

select smallFCs0, FlexClustObj and id 3087-3092+3186-3191+3365-3370+3408-3413+3426-3431+3432-3437
color brown, smallFCs0
show sticks, smallFCs0
show cartoon, smallFCs0
cartoon tube, smallFCs0



select RC, FlexClustObj and id 1-1319
color blue, RC
show lines, RC


select dangle, FlexClustObj and id 1320-3069
color white, dangle
show lines, dangle


select HPHOB3, FlexClustObj and resn "XXX" 
show spheres, HPHOB3
set sphere_scale=0.4
deselect HPHOB3
