delete RigidClustObj
delete RC*
delete smallRC*
delete FC*
delete HPHOB2
load 1bck_withH_flex_0001.pdb, RigidClustObj


select RC1, RigidClustObj and id 312-983
color red, RC1
cartoon automatic, RC1
show cartoon, RC1

select RC2, RigidClustObj and id 284-311
color cyan, RC2
cartoon automatic, RC2
show cartoon, RC2

select RC3, RigidClustObj and id 268-283
color green, RC3
cartoon automatic, RC3
show cartoon, RC3

select RC5, RigidClustObj and id 253-261
color purple, RC5
cartoon automatic, RC5
show cartoon, RC5

select smallRCs0, RigidClustObj and id 262-267+247-252+241-246+237-240+233-236+229-232+225-228+219-224+213-218+207-212+201-206+195-200+189-194+184-188+179-183+174-178+169-173+167-168+165-166+163-164+161-162+159-160+157-158+155-156+153-154+151-152+149-150+147-148+145-146+143-144+141-142+139-140+137-138+135-136+133-134+131-132+129-130+127-128+125-126+123-124+121-122+119-120+117-118+115-116+113-114
color brown, smallRCs0
show sticks, smallRCs0

select smallRCs1, RigidClustObj and id 111-112+109-110+107-108+105-106+103-104+101-102+99-100+97-98+95-96+93-94+91-92+89-90+87-88+85-86+83-84+81-82+79-80+77-78+75-76+73-74+71-72+69-70+67-68+65-66+63-64+61-62+59-60+57-58+55-56+53-54+51-52+49-50+47-48+45-46+43-44+41-42+39-40+37-38+35-36+33-34+31-32+29-30+27-28+25-26+23-24+21-22+19-20+17-18+15-16+13-14
color brown, smallRCs1
show sticks, smallRCs1

select smallRCs2, RigidClustObj and id 11-12+9-10+7-8+5-6+3-4+1-2


select FC, RigidClustObj and id 984-1948
color white, FC
show lines, FC

select HPHOB2, RigidClustObj and resn "XXX" 
show spheres, HPHOB2
set sphere_scale=0.4
deselect HPHOB2
