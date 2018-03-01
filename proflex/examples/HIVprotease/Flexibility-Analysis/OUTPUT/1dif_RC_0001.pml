delete RigidClustObj
delete RC*
delete smallRC*
delete FC*
delete HPHOB2
load 1dif_flex_0001.pdb, RigidClustObj


select RC1, RigidClustObj and id 132-1319
color red, RC1
cartoon automatic, RC1
show cartoon, RC1

select RC3, RigidClustObj and id 117-125
color green, RC3
cartoon automatic, RC3
show cartoon, RC3

select smallRCs0, RigidClustObj and id 126-131+111-116+105-110+99-104+93-98+89-92+85-88+81-84+76-80+71-75+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+71-70+69-70+69-68+67-68+65-66+63-64+61-62+59-60+57-58
color brown, smallRCs0
show sticks, smallRCs0

select smallRCs1, RigidClustObj and id 55-56+55-54+53-54+53-52+51-52+51-50+49-50+47-48+45-46+43-44+43-42+41-42+39-40+37-38+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36+37-36
color brown, smallRCs1
show sticks, smallRCs1

select smallRCs2, RigidClustObj and id 35-36+35-34+33-34+33-32+31-32+29-30+27-28+25-26+23-24+21-22+21-20+21-20+21-20+19-20+19-18+17-18+15-16+13-14+11-12+9-10+9-8+7-8+5-6+3-4+1-2


select FC, RigidClustObj and id 1320-3519
color white, FC
show lines, FC

select HPHOB2, RigidClustObj and resn "XXX" 
show spheres, HPHOB2
set sphere_scale=0.4
deselect HPHOB2
