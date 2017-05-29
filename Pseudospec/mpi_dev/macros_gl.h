#define ORDER4 1

#ifdef ORDER4
  integer, parameter :: order=2  ! order of integrator is 2*s
  real(dl), parameter :: a(order,order) = reshape( &
       (/ 0.25, 0.25-0.5/sqrt(3._dl), &
       0.25+0.5/sqrt(3._dl), 0.25/), [order,order])
  real(dl), parameter :: b(order) = (/0.5,0.5/)
#endif

!!!!! 6th order scheme !!!!!!!
#ifdef ORDER6
  integer, parameter :: order=3
  real(dl), parameter :: a(order,order) = reshape( (/ &
       5.0/36.0, 2.0/9.0-1.0/sqrt(15.0), 5.0/36.0-0.5/sqrt(15.0), &
       5.0/36.0+sqrt(15.0)/24.0, 2.0/9.0, 5.0/36.0-sqrt(15.0)/24.0, &
       5.0/36.0+0.5/sqrt(15.0), 2.0/9.0+1.0/sqrt(15.0), 5.0/36.0 /) &
       , [order, order])
  real(dl), parameter :: b(order) = (/ 5.0/18.0, 4.0/9.0, 5.0/18.0/)
#endif

!!!!! 8th order scheme !!!!!!!
#ifdef ORDER8
  integer, parameter :: order=4
  real(dl), parameter :: a(order,order) = reshape( (/ &
       0.869637112843634643432659873054998518D-1, -0.266041800849987933133851304769531093D-1, &
         0.126274626894047245150568805746180936D-1, -0.355514968579568315691098184956958860D-2, &
         0.188118117499868071650685545087171160D0,   0.163036288715636535656734012694500148D0,  &
         -0.278804286024708952241511064189974107D-1,  0.673550059453815551539866908570375889D-2, &
         0.167191921974188773171133305525295945D0,   0.353953006033743966537619131807997707D0,  &
         0.163036288715636535656734012694500148D0,  -0.141906949311411429641535704761714564D-1, &
         0.177482572254522611843442956460569292D0,   0.313445114741868346798411144814382203D0,  &
         0.352676757516271864626853155865953406D0,   0.869637112843634643432659873054998518D-1 /), [order,order] )
    real, parameter ::   b(order) = (/ &
         0.173927422568726928686531974610999704D0,   0.326072577431273071313468025389000296D0,  &
         0.326072577431273071313468025389000296D0,   0.173927422568726928686531974610999704D0  /)
#endif

!!!!! 10th order scheme !!!!!!
#ifdef ORDER10
  integer, parameter :: order=5
  real(dl), parameter :: a(order,order) = reshape( (/ &
         0.5923172126404727187856601017997934066D-1, -1.9570364359076037492643214050884060018D-2, &
         1.1254400818642955552716244215090748773D-2, -0.5593793660812184876817721964475928216D-2, &
         1.5881129678659985393652424705934162371D-3,  1.2815100567004528349616684832951382219D-1, &
         1.1965716762484161701032287870890954823D-1, -2.4592114619642200389318251686004016630D-2, &
         1.0318280670683357408953945056355839486D-2, -2.7689943987696030442826307588795957613D-3, &
         1.1377628800422460252874127381536557686D-1,  2.6000465168064151859240589518757397939D-1, &
         1.4222222222222222222222222222222222222D-1, -2.0690316430958284571760137769754882933D-2, &
         4.6871545238699412283907465445931044619D-3,  1.2123243692686414680141465111883827708D-1, &
         2.2899605457899987661169181236146325697D-1,  3.0903655906408664483376269613044846112D-1, &
         1.1965716762484161701032287870890954823D-1, -0.9687563141950739739034827969555140871D-2, &
         1.1687532956022854521776677788936526508D-1,  2.4490812891049541889746347938229502468D-1, &
         2.7319004362580148889172820022935369566D-1,  2.5888469960875927151328897146870315648D-1, &
         0.5923172126404727187856601017997934066D-1 /) , [order,order])
    real, parameter :: b(order) = (/ &
         1.1846344252809454375713202035995868132D-1,  2.3931433524968323402064575741781909646D-1, &
         2.8444444444444444444444444444444444444D-1,  2.3931433524968323402064575741781909646D-1, &
         1.1846344252809454375713202035995868132D-1 /)
#endif