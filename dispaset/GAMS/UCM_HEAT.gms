$ontext
This file is coupling the UCM with the DSM heat model by the mean of a state
and restart statement.
$offtext

$set gamsparm "ide=%gams.ide% lo=%gams.lo% errorlog=%gams.errorlog% errmsg=1"
$show

execute "gams DSM.gms %gamsparm% s=DSM"
execute "gams UCM_h.gms %gamsparm% r=DSM s=UCM"

