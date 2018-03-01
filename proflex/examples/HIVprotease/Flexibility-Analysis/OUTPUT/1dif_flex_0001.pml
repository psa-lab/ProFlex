delete FlexIndexObj
delete HPHOB1
load 1dif_flex_0001.pdb, FlexIndexObj

cmd.color( 'blue', 'FlexIndexObj and b< 42.0')
cmd.color('tv_blue', 'FlexIndexObj and b> 42.0 and b < 44.0')
cmd.color('tv_blue', 'FlexIndexObj and b= 44.0')
cmd.color('marine', 'FlexIndexObj and b> 44.0 and b < 46.0')
cmd.color('marine', 'FlexIndexObj and b= 46.0')
cmd.color('lightblue', 'FlexIndexObj and b> 46.0 and b < 48.0')
cmd.color('lightblue', 'FlexIndexObj and b= 48.0')
cmd.color('cyan', 'FlexIndexObj and b> 48.0 and b < 49.0')
cmd.color('cyan', 'FlexIndexObj and b= 49.0')
cmd.color('grey', 'FlexIndexObj and b> 49.0 and b < 50.0')
cmd.color('grey', 'FlexIndexObj and b= 50.0')
cmd.color('yellow', 'FlexIndexObj and b> 50.0 and b < 52.0')
cmd.color('yellow', 'FlexIndexObj and b= 52.0')
cmd.color('yelloworange', 'FlexIndexObj and b> 52.0 and b < 54.0')
cmd.color('yelloworange', 'FlexIndexObj and b= 54.0')
cmd.color('orange', 'FlexIndexObj and b> 54.0 and b < 56.0')
cmd.color('orange', 'FlexIndexObj and b= 56.0')
cmd.color('tv_red', 'FlexIndexObj and b> 56.0 and b < 58.0')
cmd.color('tv_red', 'FlexIndexObj and b= 58.0')
cmd.color('red', 'FlexIndexObj and b> 58')


select HPHOB1, FlexIndexObj and resn "XXX" 
show spheres, HPHOB1
set sphere_scale=0.4
deselect HPHOB1
