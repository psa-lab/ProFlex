delete FlexClustObj
delete FC*
delete smallFC*
delete RC*
delete HPHOB3
load 1bck_withH_flex_0001.pdb, FlexClustObj


select FC1, FlexClustObj and id 1499-1936
color wheat, FC1
show sticks, FC1
show cartoon, FC1
cartoon tube, FC1

select smallFCs0, FlexClustObj and id 1937-1942+1943-1948
color brown, smallFCs0
show sticks, smallFCs0
show cartoon, smallFCs0
cartoon tube, smallFCs0



select RC, FlexClustObj and id 1-983
color blue, RC
show lines, RC


select dangle, FlexClustObj and id 984-1498
color white, dangle
show lines, dangle


select HPHOB3, FlexClustObj and resn "XXX" 
show spheres, HPHOB3
set sphere_scale=0.4
deselect HPHOB3
