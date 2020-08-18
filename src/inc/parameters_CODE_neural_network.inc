C======================================================================|
C----------------------------------------------------------------------|
C     Parameters for the Dreicer neural network
C----------------------------------------------------------------------|

      real(kind=dp), parameter, dimension(8) ::
     >  input_std = (/ 
     >   4.9612189726200_dp,   5.5830833950300_dp,   0.1953846218470_dp,
     >   4.9828199651400_dp,   2.4518547235100_dp,   0.3040905264250_dp,
     >   0.0345813696948_dp,   2.2886290216900_dp
     >      /)

      real(kind=dp), parameter, dimension(8) ::
     >  input_mean = (/ 
     >   7.1561340244000_dp,   9.2037474934600_dp,   0.1059798956440_dp,
     >   7.0891001224100_dp,  47.8203654981000_dp,   0.5363477883270_dp,
     >   0.0699809427188_dp,  -9.4087288982000_dp
     >      /)

      real(kind=dp), parameter ::
     >  output_std =   18.5581667797000_dp

      real(kind=dp), parameter ::
     >  output_mean =  -23.1854797433000_dp

      real(kind=dp), parameter, dimension(20) ::
     >  b1 = (/ 
     >   0.9228090000000_dp,  -0.6333320000000_dp,   0.9568790000000_dp,
     >   0.5677170000000_dp,   1.1156800000000_dp,   0.6332590000000_dp,
     >   0.7142460000000_dp,   0.4260980000000_dp,   0.5753610000000_dp,
     >   0.5778780000000_dp,   0.9363330000000_dp,  -0.6446000000000_dp,
     >  -0.4482020000000_dp,   0.5006270000000_dp,  -0.3222090000000_dp,
     >  -0.3031710000000_dp,  -0.6810620000000_dp,  -0.2193820000000_dp,
     >   0.3405370000000_dp,   0.5261250000000_dp
     >      /)

      real(kind=dp), parameter, dimension(20,8) ::
     >  W1 = reshape( (/ 
     >  -0.0300745000000_dp,  -0.0057330200000_dp,  -0.0071478800000_dp,
     >   0.0445892000000_dp,  -0.0759706000000_dp,  -0.1283420000000_dp,
     >   0.0269302000000_dp,  -0.2294620000000_dp,  -0.1026160000000_dp,
     >  -0.1140890000000_dp,   0.1580920000000_dp,   0.0903593000000_dp,
     >   0.2207100000000_dp,  -0.2398530000000_dp,  -0.1236780000000_dp,
     >   0.5379190000000_dp,   0.5227460000000_dp,  -0.4120660000000_dp,
     >   0.5940480000000_dp,  -0.0101526000000_dp,  -0.0492966000000_dp,
     >   0.1087650000000_dp,   0.0555219000000_dp,  -0.0060380000000_dp,
     >   0.1338720000000_dp,   0.0882125000000_dp,  -0.1475170000000_dp,
     >   0.3388990000000_dp,   0.0592563000000_dp,  -0.2334690000000_dp,
     >  -0.1120470000000_dp,  -0.0002139150000_dp,   0.0827269000000_dp,
     >   0.1288950000000_dp,  -0.2306690000000_dp,   0.2639790000000_dp,
     >   0.0748574000000_dp,   0.1376710000000_dp,  -0.1406740000000_dp,
     >   0.0160230000000_dp,  -0.0259351000000_dp,  -0.0796380000000_dp,
     >  -0.0839504000000_dp,   0.0154224000000_dp,   0.6576040000000_dp,
     >  -0.1230200000000_dp,  -0.1227260000000_dp,   0.7406510000000_dp,
     >   0.2322500000000_dp,  -0.0731540000000_dp,  -0.0659573000000_dp,
     >  -0.0953753000000_dp,   0.1480280000000_dp,   0.0402638000000_dp,
     >   0.3899360000000_dp,  -0.0048738600000_dp,  -0.5980240000000_dp,
     >  -0.0102158000000_dp,   0.1029230000000_dp,  -0.0251902000000_dp,
     >  -0.1977780000000_dp,  -0.0446178000000_dp,   0.1662420000000_dp,
     >   0.1874880000000_dp,  -0.0339612000000_dp,   0.0638418000000_dp,
     >  -0.4950900000000_dp,  -0.1995370000000_dp,   0.4623550000000_dp,
     >   0.2049250000000_dp,   0.0238271000000_dp,  -0.0851210000000_dp,
     >  -0.0947850000000_dp,   0.2416350000000_dp,   0.2216660000000_dp,
     >  -0.2941190000000_dp,  -0.2940860000000_dp,   0.5702880000000_dp,
     >  -0.0371824000000_dp,   0.0382110000000_dp,   0.0074670200000_dp,
     >   0.0284712000000_dp,   0.0119449000000_dp,  -0.0184206000000_dp,
     >   0.0072863900000_dp,  -0.2200960000000_dp,   0.0044491000000_dp,
     >   0.0154031000000_dp,  -0.0017818000000_dp,   0.0608162000000_dp,
     >  -0.0157809000000_dp,   0.1236120000000_dp,  -0.0225893000000_dp,
     >   0.0516254000000_dp,  -0.0983403000000_dp,  -0.0167992000000_dp,
     >   0.0027014900000_dp,   0.0135675000000_dp,  -0.0358325000000_dp,
     >  -0.0152584000000_dp,  -0.2924820000000_dp,   0.3698770000000_dp,
     >   0.0110446000000_dp,   0.3103830000000_dp,   0.6346320000000_dp,
     >   0.2291340000000_dp,  -0.4549550000000_dp,   0.3344110000000_dp,
     >   0.3639350000000_dp,  -0.2442640000000_dp,   0.5497090000000_dp,
     >  -0.3408240000000_dp,  -0.6835860000000_dp,   0.3562470000000_dp,
     >  -0.2519970000000_dp,   0.0143781000000_dp,  -0.5078660000000_dp,
     >  -0.3180620000000_dp,  -0.3624270000000_dp,   0.2110740000000_dp,
     >   0.4864790000000_dp,  -0.5736680000000_dp,  -0.4182850000000_dp,
     >   0.1910480000000_dp,  -0.0816593000000_dp,   0.0219170000000_dp,
     >   0.5301440000000_dp,  -0.0424001000000_dp,  -0.0187410000000_dp,
     >  -0.1740560000000_dp,  -0.3560530000000_dp,   0.1977580000000_dp,
     >  -0.1450010000000_dp,   0.4149200000000_dp,   0.1759460000000_dp,
     >   0.2136220000000_dp,  -0.0245604000000_dp,   0.0301480000000_dp,
     >   0.0683319000000_dp,   0.6641830000000_dp,  -0.2117720000000_dp,
     >  -0.2766340000000_dp,  -0.4056130000000_dp,   0.4121880000000_dp,
     >   0.1268240000000_dp,  -0.1573440000000_dp,   0.4463650000000_dp,
     >   0.4110410000000_dp,  -0.6920860000000_dp,   0.1539580000000_dp,
     >   0.0474924000000_dp,  -0.3962310000000_dp,   0.5079440000000_dp,
     >  -0.3784760000000_dp,  -0.4079280000000_dp,  -0.1887740000000_dp,
     >   0.1659610000000_dp,   0.3542690000000_dp,  -0.1010850000000_dp,
     >   0.5266140000000_dp
     >      /), (/20,8/))

      real(kind=dp), parameter, dimension(20) ::
     >  b2 = (/ 
     >   0.1737450000000_dp,   0.6120800000000_dp,   0.3439650000000_dp,
     >  -0.1308900000000_dp,   0.5176530000000_dp,   0.7533340000000_dp,
     >   0.3806170000000_dp,  -0.2164450000000_dp,   0.5006960000000_dp,
     >  -0.1386610000000_dp,   0.3067670000000_dp,  -0.7422290000000_dp,
     >  -0.8370640000000_dp,   0.2174200000000_dp,   0.2782600000000_dp,
     >  -0.6934510000000_dp,   0.5551420000000_dp,   0.4280250000000_dp,
     >  -0.0885804000000_dp,   0.0908920000000_dp
     >      /)

      real(kind=dp), parameter, dimension(20,20) ::
     >  W2 = reshape( (/ 
     >  -0.5863500000000_dp,  -0.4225920000000_dp,  -0.1025900000000_dp,
     >   0.1911370000000_dp,  -0.1171500000000_dp,   0.2470980000000_dp,
     >  -0.5644890000000_dp,   0.3555410000000_dp,  -0.6810300000000_dp,
     >   0.7109330000000_dp,  -0.3101680000000_dp,   0.5402860000000_dp,
     >   0.2234270000000_dp,   0.3333960000000_dp,  -0.6137060000000_dp,
     >   0.3022660000000_dp,  -0.3268770000000_dp,   0.1745680000000_dp,
     >   0.8160780000000_dp,  -0.7015100000000_dp,  -0.4386720000000_dp,
     >   0.0521457000000_dp,   0.2978710000000_dp,  -0.2680880000000_dp,
     >   0.3695240000000_dp,  -0.0310993000000_dp,   0.5223820000000_dp,
     >   0.4660000000000_dp,  -0.0243779000000_dp,  -0.2139240000000_dp,
     >  -0.1369760000000_dp,   0.2160310000000_dp,  -0.4905140000000_dp,
     >   0.3003190000000_dp,   0.4961610000000_dp,  -0.0420151000000_dp,
     >  -0.0859462000000_dp,  -0.1494240000000_dp,  -0.1143460000000_dp,
     >   0.0997946000000_dp,  -0.6038800000000_dp,   0.0376996000000_dp,
     >   0.2969500000000_dp,   0.5364500000000_dp,  -0.1837790000000_dp,
     >  -0.0253376000000_dp,   0.1383140000000_dp,   0.5910470000000_dp,
     >  -0.1798500000000_dp,  -0.4466440000000_dp,  -0.1153310000000_dp,
     >   0.4514830000000_dp,  -0.5661730000000_dp,   0.8273960000000_dp,
     >  -0.3693210000000_dp,  -0.3382040000000_dp,  -0.3858100000000_dp,
     >   0.2337080000000_dp,   0.2187930000000_dp,  -0.3854320000000_dp,
     >  -0.5614890000000_dp,  -0.2356700000000_dp,  -0.0008624720000_dp,
     >  -0.0632909000000_dp,  -0.1014540000000_dp,   0.7627470000000_dp,
     >   0.2700940000000_dp,   0.0420865000000_dp,  -0.4266740000000_dp,
     >   0.1812180000000_dp,   0.2808810000000_dp,   0.5059940000000_dp,
     >  -0.1415920000000_dp,   0.2655020000000_dp,  -0.3622020000000_dp,
     >   0.7819860000000_dp,  -0.4408990000000_dp,   0.5915140000000_dp,
     >   0.1688120000000_dp,   0.3558370000000_dp,  -0.3644130000000_dp,
     >   0.0030770300000_dp,  -0.0799795000000_dp,  -0.4615560000000_dp,
     >   0.0506903000000_dp,   1.6756500000000_dp,  -0.1108040000000_dp,
     >   0.0617810000000_dp,  -0.2988010000000_dp,   0.2414630000000_dp,
     >   2.8216800000000_dp,  -0.2890620000000_dp,  -0.0121970000000_dp,
     >   0.6769450000000_dp,   0.2828920000000_dp,   0.0558212000000_dp,
     >  -0.1005170000000_dp,   1.4483100000000_dp,   0.0030488000000_dp,
     >  -0.0665813000000_dp,   0.0203549000000_dp,   0.1508210000000_dp,
     >  -0.1404290000000_dp,   0.3317220000000_dp,  -0.1866280000000_dp,
     >   0.2434750000000_dp,   0.1182430000000_dp,   0.2211420000000_dp,
     >   0.1083370000000_dp,  -0.3385750000000_dp,   0.2735550000000_dp,
     >  -0.1793540000000_dp,  -0.4322670000000_dp,   0.4457130000000_dp,
     >  -0.2236040000000_dp,   0.1743490000000_dp,   0.1559860000000_dp,
     >   0.2994330000000_dp,   0.3114660000000_dp,  -0.0892183000000_dp,
     >  -0.1595760000000_dp,  -0.3954770000000_dp,  -0.4513220000000_dp,
     >   0.1742810000000_dp,  -0.2103270000000_dp,   0.0943831000000_dp,
     >   0.0030608900000_dp,  -0.0226697000000_dp,  -0.2966740000000_dp,
     >  -0.4749370000000_dp,  -0.3794020000000_dp,  -0.1666580000000_dp,
     >   0.2116630000000_dp,  -0.4305020000000_dp,  -0.1718220000000_dp,
     >  -0.1841150000000_dp,  -0.2981940000000_dp,   0.1039490000000_dp,
     >   0.2724320000000_dp,   0.2673900000000_dp,  -0.1122070000000_dp,
     >   0.1577410000000_dp,   0.2186210000000_dp,  -0.2218750000000_dp,
     >   0.0071987300000_dp,   1.0586400000000_dp,  -0.3470420000000_dp,
     >  -0.0946886000000_dp,  -0.2384580000000_dp,  -0.3167320000000_dp,
     >   1.3799300000000_dp,  -0.0761549000000_dp,  -0.5029160000000_dp,
     >  -0.1742290000000_dp,  -0.1213660000000_dp,  -0.1411580000000_dp,
     >   0.1081260000000_dp,  -0.0791139000000_dp,   0.2346850000000_dp,
     >  -0.1509910000000_dp,  -0.3403000000000_dp,  -0.1582230000000_dp,
     >  -0.1104910000000_dp,  -0.0629570000000_dp,  -0.6405360000000_dp,
     >   0.4763290000000_dp,   0.2456920000000_dp,   0.5940800000000_dp,
     >   0.1720140000000_dp,  -0.1151340000000_dp,  -0.0863799000000_dp,
     >   0.0110979000000_dp,   0.2493830000000_dp,   0.6757190000000_dp,
     >  -0.5602290000000_dp,   0.0937934000000_dp,  -0.4382610000000_dp,
     >   0.5024670000000_dp,   0.1312890000000_dp,   0.2001240000000_dp,
     >   0.0022411200000_dp,   0.0237450000000_dp,   0.0777079000000_dp,
     >  -0.0766876000000_dp,   0.2937780000000_dp,   0.2316290000000_dp,
     >   0.2864230000000_dp,  -0.2101480000000_dp,  -0.1748080000000_dp,
     >   0.1088250000000_dp,   0.1523670000000_dp,  -0.0506017000000_dp,
     >  -0.5807260000000_dp,   0.1463020000000_dp,   0.1910050000000_dp,
     >   0.0356107000000_dp,  -0.0530153000000_dp,   0.1992240000000_dp,
     >  -0.1159500000000_dp,  -0.0917435000000_dp,  -0.4532410000000_dp,
     >   0.0592980000000_dp,  -0.3548040000000_dp,  -0.0595797000000_dp,
     >  -0.1235770000000_dp,   0.5712260000000_dp,   0.2978220000000_dp,
     >   0.7634980000000_dp,  -0.4603590000000_dp,   0.0994124000000_dp,
     >  -0.0149924000000_dp,  -0.1413910000000_dp,   0.3296310000000_dp,
     >   0.6146890000000_dp,   0.2694530000000_dp,  -0.3575700000000_dp,
     >   0.1587980000000_dp,   0.6662860000000_dp,   0.0434260000000_dp,
     >   0.2610950000000_dp,  -0.1031280000000_dp,  -0.0128884000000_dp,
     >   0.1943330000000_dp,   0.4388840000000_dp,  -0.2941180000000_dp,
     >  -0.3634620000000_dp,  -0.1932180000000_dp,  -0.1343220000000_dp,
     >   0.2592880000000_dp,  -0.4306610000000_dp,  -0.0248687000000_dp,
     >  -0.0692775000000_dp,   0.0119444000000_dp,  -0.3015440000000_dp,
     >   0.3427750000000_dp,  -0.0977751000000_dp,  -0.1023000000000_dp,
     >  -0.2808470000000_dp,   0.0842895000000_dp,   0.1502450000000_dp,
     >  -0.3560280000000_dp,   0.0485015000000_dp,   0.4631750000000_dp,
     >  -0.1919200000000_dp,   0.6237910000000_dp,  -0.6004020000000_dp,
     >   0.1444550000000_dp,  -0.4218660000000_dp,  -0.0853164000000_dp,
     >  -0.3984330000000_dp,  -0.2180920000000_dp,   0.3351160000000_dp,
     >  -0.2851950000000_dp,  -1.3306100000000_dp,   0.0880120000000_dp,
     >  -0.2857800000000_dp,   0.0648254000000_dp,  -0.2404660000000_dp,
     >  -0.0603468000000_dp,  -0.3323180000000_dp,   0.0308356000000_dp,
     >   0.2631340000000_dp,  -0.3748280000000_dp,   0.2382880000000_dp,
     >  -0.3018240000000_dp,   0.5441630000000_dp,   0.7320840000000_dp,
     >  -0.1451220000000_dp,  -0.2585290000000_dp,   0.0478932000000_dp,
     >   0.5198530000000_dp,  -0.1212170000000_dp,  -0.0966694000000_dp,
     >   0.1803190000000_dp,  -0.3975680000000_dp,   0.4028870000000_dp,
     >   0.0871810000000_dp,   0.3362040000000_dp,   0.1435230000000_dp,
     >  -0.5377080000000_dp,   0.0334945000000_dp,   0.1566030000000_dp,
     >   0.0936719000000_dp,  -0.6123820000000_dp,  -0.1486300000000_dp,
     >   0.1647790000000_dp,   0.0165015000000_dp,   0.4382180000000_dp,
     >   0.0985250000000_dp,  -0.2768520000000_dp,   0.7558260000000_dp,
     >  -0.0360731000000_dp,   0.6917130000000_dp,  -0.5064710000000_dp,
     >  -0.3100830000000_dp,   0.3738950000000_dp,  -0.0706768000000_dp,
     >  -0.2330240000000_dp,  -0.0036990200000_dp,   0.2632680000000_dp,
     >  -0.1947220000000_dp,  -0.1958110000000_dp,  -0.1779570000000_dp,
     >   0.4590690000000_dp,  -0.3666140000000_dp,  -0.5638330000000_dp,
     >   0.3761140000000_dp,  -0.0078393200000_dp,   0.0729553000000_dp,
     >   0.1886540000000_dp,  -1.4147600000000_dp,  -0.1656310000000_dp,
     >  -0.1680060000000_dp,   0.0911605000000_dp,   0.0097567500000_dp,
     >   0.0522160000000_dp,   0.0085970800000_dp,  -0.3969380000000_dp,
     >   0.3147390000000_dp,   0.1064770000000_dp,  -0.0678543000000_dp,
     >   0.1529930000000_dp,   0.0195798000000_dp,  -0.3079960000000_dp,
     >   0.2113380000000_dp,  -3.5275700000000_dp,  -0.2241390000000_dp,
     >  -0.2313990000000_dp,  -0.0739903000000_dp,  -0.3683100000000_dp,
     >  -3.1158300000000_dp,   0.0378446000000_dp,   0.0503492000000_dp,
     >  -0.4184290000000_dp,  -0.0031181000000_dp,   0.1672670000000_dp,
     >  -0.1301950000000_dp,  -2.6622500000000_dp,   0.2609850000000_dp,
     >  -0.1013420000000_dp,   0.4872540000000_dp,   0.2686130000000_dp,
     >   0.2025700000000_dp,  -0.7042760000000_dp,  -0.4666370000000_dp,
     >   0.0706474000000_dp,   0.1091170000000_dp,   0.0323285000000_dp,
     >  -0.3786660000000_dp,  -0.2512000000000_dp,  -0.1055130000000_dp,
     >   0.0369415000000_dp,   0.1842510000000_dp,   0.0248688000000_dp,
     >   0.0717479000000_dp,  -0.0270453000000_dp,  -0.2692260000000_dp,
     >  -0.2711120000000_dp,   0.1054070000000_dp,  -0.1082450000000_dp,
     >   0.2333830000000_dp,   0.0971745000000_dp,  -0.1660480000000_dp,
     >  -0.1590860000000_dp,  -0.3451650000000_dp,  -0.2836590000000_dp,
     >   0.1084500000000_dp,  -0.3165760000000_dp,  -0.0755463000000_dp,
     >   0.1975300000000_dp,   0.2927120000000_dp,   0.1858390000000_dp,
     >  -0.4565300000000_dp,   0.3267930000000_dp,  -0.0646299000000_dp,
     >  -0.3615320000000_dp,   0.0466244000000_dp,   0.0651703000000_dp,
     >   0.1598630000000_dp,   0.3276840000000_dp,   0.7287130000000_dp,
     >  -0.2037810000000_dp,  -0.5255940000000_dp,   0.2953160000000_dp,
     >   0.1978770000000_dp,   0.3092830000000_dp,   0.2254800000000_dp,
     >   0.3413610000000_dp,  -0.0974758000000_dp,   0.6037210000000_dp,
     >  -0.1047420000000_dp,   0.3296590000000_dp,  -0.2224120000000_dp,
     >  -0.0338085000000_dp,  -0.2090460000000_dp,   0.0365815000000_dp,
     >  -0.2624940000000_dp,  -0.0341975000000_dp,   0.6698640000000_dp,
     >  -0.5571630000000_dp
     >      /), (/20,20/))

      real(kind=dp), parameter, dimension(20) ::
     >  b3 = (/ 
     >   0.1281720000000_dp,  -0.6836510000000_dp,  -0.5169390000000_dp,
     >   0.7455950000000_dp,  -0.4886220000000_dp,   0.1317110000000_dp,
     >  -0.5170390000000_dp,   0.4694130000000_dp,  -0.1995410000000_dp,
     >  -0.5877640000000_dp,  -0.6841570000000_dp,  -0.1824670000000_dp,
     >  -0.3704320000000_dp,   0.2440920000000_dp,  -0.5825570000000_dp,
     >   0.7021700000000_dp,  -0.5111290000000_dp,   0.5139630000000_dp,
     >   0.3058890000000_dp,   0.5028790000000_dp
     >      /)

      real(kind=dp), parameter, dimension(20,20) ::
     >  W3 = reshape( (/ 
     >   0.7986480000000_dp,  -0.2713280000000_dp,  -0.5151500000000_dp,
     >  -0.2096990000000_dp,  -0.1627080000000_dp,   0.4871920000000_dp,
     >   0.2879370000000_dp,   0.2015580000000_dp,   0.3942950000000_dp,
     >  -0.7538470000000_dp,  -0.7315390000000_dp,  -0.1191070000000_dp,
     >   0.7844170000000_dp,  -0.1192050000000_dp,  -0.1789620000000_dp,
     >   1.0819800000000_dp,  -0.6693740000000_dp,   0.0364462000000_dp,
     >   1.1495900000000_dp,   0.3426660000000_dp,  -0.5902860000000_dp,
     >  -0.1425950000000_dp,  -0.3116380000000_dp,   0.6937540000000_dp,
     >   0.3728700000000_dp,   0.3646390000000_dp,  -0.0990251000000_dp,
     >  -0.4132540000000_dp,  -0.1215470000000_dp,   0.0625514000000_dp,
     >  -0.6343460000000_dp,  -0.8453740000000_dp,  -0.0039153500000_dp,
     >   0.0378652000000_dp,  -0.7239670000000_dp,   0.6435950000000_dp,
     >  -0.3818540000000_dp,   0.3259020000000_dp,   0.1750750000000_dp,
     >   0.5878930000000_dp,   0.1616050000000_dp,  -0.2993670000000_dp,
     >  -0.3979010000000_dp,   0.7807320000000_dp,  -0.2542900000000_dp,
     >   0.1616520000000_dp,  -0.5753130000000_dp,   0.3022510000000_dp,
     >  -0.6944810000000_dp,   0.8553960000000_dp,  -0.3748710000000_dp,
     >  -0.7585080000000_dp,  -0.4500460000000_dp,   0.6623980000000_dp,
     >  -0.3469120000000_dp,   1.7273500000000_dp,   0.4415820000000_dp,
     >   0.3991080000000_dp,  -0.0465164000000_dp,   1.1002200000000_dp,
     >  -0.2225020000000_dp,   0.1615580000000_dp,   1.3447400000000_dp,
     >  -0.0140950000000_dp,   0.1192340000000_dp,  -0.1617190000000_dp,
     >   0.0259689000000_dp,  -0.6029310000000_dp,  -0.6249350000000_dp,
     >   0.6414910000000_dp,  -0.0128045000000_dp,  -0.0598682000000_dp,
     >  -0.5424490000000_dp,  -0.1551410000000_dp,   0.3727620000000_dp,
     >  -0.1813530000000_dp,   0.1359910000000_dp,  -0.2518990000000_dp,
     >  -0.5337360000000_dp,  -0.2320900000000_dp,   1.4879200000000_dp,
     >   0.0620801000000_dp,  -0.4404520000000_dp,   0.1493670000000_dp,
     >  -1.5145800000000_dp,   0.8710220000000_dp,   0.0172405000000_dp,
     >   0.5173950000000_dp,  -0.2724480000000_dp,   0.3492130000000_dp,
     >  -0.5867630000000_dp,   0.0835409000000_dp,   0.0777790000000_dp,
     >   0.3371160000000_dp,   0.4295240000000_dp,   1.9812500000000_dp,
     >   0.3519110000000_dp,   0.3813530000000_dp,  -0.1105470000000_dp,
     >   0.3959670000000_dp,  -5.1758300000000_dp,   0.2446120000000_dp,
     >   1.6142700000000_dp,   0.9874560000000_dp,   1.0017000000000_dp,
     >  -8.2865600000000_dp,   0.3810920000000_dp,   0.2708020000000_dp,
     >  -0.1120720000000_dp,   1.1456100000000_dp,   0.7131750000000_dp,
     >   0.3563160000000_dp,   0.3715690000000_dp,  -0.0292891000000_dp,
     >  -0.6942080000000_dp,  -0.6145220000000_dp,   0.2213030000000_dp,
     >  -0.8197040000000_dp,  -0.9533430000000_dp,  -0.7171700000000_dp,
     >  -0.7821660000000_dp,  -0.1703880000000_dp,  -0.2791030000000_dp,
     >  -0.0496421000000_dp,   0.7190260000000_dp,   2.0282700000000_dp,
     >  -0.3323100000000_dp,   0.7710040000000_dp,   0.5515070000000_dp,
     >  -0.2798420000000_dp,   0.0406906000000_dp,   0.7313140000000_dp,
     >  -0.1073450000000_dp,  -0.6098130000000_dp,  -1.8616000000000_dp,
     >  -0.8115150000000_dp,  -0.0537429000000_dp,  -0.4022310000000_dp,
     >   0.2226380000000_dp,   0.1803470000000_dp,  -0.0362775000000_dp,
     >   0.1535380000000_dp,   0.8965920000000_dp,  -0.0203096000000_dp,
     >   0.4081860000000_dp,   0.8209190000000_dp,   0.0615488000000_dp,
     >   0.0933618000000_dp,   0.0401361000000_dp,   0.6051780000000_dp,
     >   0.9132330000000_dp,  -0.0835606000000_dp,   0.0089849900000_dp,
     >  -0.1301050000000_dp,   0.8078060000000_dp,  -0.9396710000000_dp,
     >   0.5473740000000_dp,  -0.3250220000000_dp,  -0.4540970000000_dp,
     >  -0.8123770000000_dp,  -0.2952840000000_dp,  -0.4575580000000_dp,
     >  -1.5159900000000_dp,   0.1529740000000_dp,   0.7888350000000_dp,
     >   0.8921920000000_dp,  -0.1696350000000_dp,  -0.2462830000000_dp,
     >  -0.0859958000000_dp,  -1.4591700000000_dp,  -1.2792200000000_dp,
     >  -1.1030800000000_dp,  -0.0925234000000_dp,  -0.1747040000000_dp,
     >  -0.2630040000000_dp,   1.3227900000000_dp,  -0.9528250000000_dp,
     >   0.9042040000000_dp,   0.6533190000000_dp,   0.5226630000000_dp,
     >  -0.6311990000000_dp,   0.5335220000000_dp,   0.7718180000000_dp,
     >   0.3290640000000_dp,  -0.7992410000000_dp,  -0.9852000000000_dp,
     >  -0.2064180000000_dp,   0.3860880000000_dp,   0.2065280000000_dp,
     >   0.6255700000000_dp,   0.1786580000000_dp,  -0.5147210000000_dp,
     >  -0.2754880000000_dp,   0.9914140000000_dp,  -0.0202029000000_dp,
     >   0.1162950000000_dp,   0.8190350000000_dp,  -0.7665210000000_dp,
     >  -0.9361250000000_dp,  -0.3402930000000_dp,  -0.7097890000000_dp,
     >   0.0670842000000_dp,  -1.0247900000000_dp,   0.0436989000000_dp,
     >  -0.2447300000000_dp,  -0.8364360000000_dp,   0.0095278400000_dp,
     >  -0.1244790000000_dp,  -0.0832871000000_dp,  -0.0161334000000_dp,
     >  -0.0311367000000_dp,  -0.0188483000000_dp,   0.3965160000000_dp,
     >   0.1282440000000_dp,  -0.5894620000000_dp,   0.3210260000000_dp,
     >  -0.2472360000000_dp,   0.1577630000000_dp,   0.3992310000000_dp,
     >   0.5136230000000_dp,   0.0869203000000_dp,   0.3686090000000_dp,
     >   1.6351100000000_dp,  -0.3612190000000_dp,  -0.6644030000000_dp,
     >   1.1059200000000_dp,  -0.4951520000000_dp,   0.1442670000000_dp,
     >   0.0623015000000_dp,   2.1122600000000_dp,   1.2926400000000_dp,
     >   0.1928870000000_dp,  -0.2702900000000_dp,   0.4931080000000_dp,
     >  -0.6495370000000_dp,  -0.9105520000000_dp,   1.3406800000000_dp,
     >  -1.2086100000000_dp,  -1.2968000000000_dp,  -0.8622960000000_dp,
     >  -0.4871620000000_dp,   0.6111920000000_dp,   0.8977540000000_dp,
     >   0.0474346000000_dp,   1.3270100000000_dp,  -0.5722460000000_dp,
     >   0.3361050000000_dp,  -0.2845180000000_dp,   0.3521370000000_dp,
     >   0.2040220000000_dp,   0.7793800000000_dp,  -0.2673960000000_dp,
     >  -0.0029137900000_dp,  -0.5875360000000_dp,   0.3951330000000_dp,
     >  -1.1477800000000_dp,   0.5376770000000_dp,  -0.6016390000000_dp,
     >  -0.0781300000000_dp,  -0.8739150000000_dp,  -0.2649810000000_dp,
     >  -0.3926840000000_dp,   0.8127510000000_dp,  -0.1147370000000_dp,
     >   0.8939460000000_dp,  -1.1588900000000_dp,   0.3093710000000_dp,
     >   0.1281030000000_dp,   0.3925310000000_dp,   0.5356000000000_dp,
     >   0.9461260000000_dp,   0.4442040000000_dp,  -0.3259250000000_dp,
     >  -0.4970310000000_dp,   0.5819770000000_dp,  -2.0608300000000_dp,
     >   0.4646710000000_dp,  -0.4075240000000_dp,  -0.8544010000000_dp,
     >  -0.8377090000000_dp,   1.6445500000000_dp,  -0.4816350000000_dp,
     >  -1.1629500000000_dp,   0.6351920000000_dp,   0.0731557000000_dp,
     >   2.2559300000000_dp,  -0.2500350000000_dp,  -0.1258920000000_dp,
     >  -0.0748858000000_dp,  -0.9062710000000_dp,  -0.3832310000000_dp,
     >   0.5368970000000_dp,   0.5440390000000_dp,   0.0968684000000_dp,
     >   0.4708280000000_dp,   0.9601210000000_dp,  -0.9286780000000_dp,
     >   0.6969190000000_dp,   0.2433190000000_dp,   1.0925000000000_dp,
     >  -0.4029630000000_dp,   0.2251430000000_dp,   1.0147800000000_dp,
     >  -0.2322970000000_dp,   0.1498160000000_dp,  -0.8642460000000_dp,
     >   0.2617390000000_dp,  -0.2911060000000_dp,   0.1774830000000_dp,
     >   0.9333250000000_dp,   0.5708160000000_dp,  -0.5307310000000_dp,
     >  -0.4330250000000_dp,   0.2680730000000_dp,  -0.1514830000000_dp,
     >  -1.0118600000000_dp,   0.0336602000000_dp,  -0.5357810000000_dp,
     >  -0.4996740000000_dp,  -0.3224010000000_dp,  -0.8892360000000_dp,
     >   0.2213310000000_dp,  -0.3033550000000_dp,  -0.1516190000000_dp,
     >  -1.3095300000000_dp,  -1.0407100000000_dp,   0.0043938600000_dp,
     >  -0.7699920000000_dp,   0.3272350000000_dp,  -0.6244370000000_dp,
     >  -0.6942620000000_dp,  -0.0020533900000_dp,   0.1173560000000_dp,
     >   0.1642210000000_dp,   0.2718750000000_dp,   1.5958500000000_dp,
     >  -0.0418462000000_dp,  -0.4037880000000_dp,   0.7281090000000_dp,
     >   0.2923710000000_dp,  -1.3105400000000_dp,  -0.0108845000000_dp,
     >   1.1142300000000_dp,  -0.0548446000000_dp,   0.0420132000000_dp,
     >  -2.2661600000000_dp,   0.0414084000000_dp,   0.0520541000000_dp,
     >   0.1152210000000_dp,   0.3707890000000_dp,   0.3467540000000_dp,
     >  -0.0394782000000_dp,  -0.3890110000000_dp,   0.1177630000000_dp,
     >   0.2310350000000_dp,  -1.2578400000000_dp,   0.0085766500000_dp,
     >  -0.3532890000000_dp,  -0.3849100000000_dp,  -0.5138590000000_dp,
     >  -0.1173740000000_dp,  -0.2005870000000_dp,   0.6623330000000_dp,
     >  -0.1241100000000_dp,   0.5758930000000_dp,  -0.0031488800000_dp,
     >  -0.1278340000000_dp,  -0.4656360000000_dp,   0.2590110000000_dp,
     >   0.4644060000000_dp,   0.6064610000000_dp,   0.3358590000000_dp,
     >  -0.3668780000000_dp,  -0.0338622000000_dp,  -0.0284809000000_dp,
     >  -1.6133300000000_dp,   0.2766290000000_dp,  -0.3684830000000_dp,
     >  -0.1473760000000_dp,  -0.1314780000000_dp,  -0.1215900000000_dp,
     >  -0.0354804000000_dp,  -0.3782560000000_dp,   0.4615860000000_dp,
     >  -0.0605411000000_dp,   0.6309250000000_dp,  -0.2998520000000_dp,
     >   0.1792970000000_dp,  -0.4958220000000_dp,  -0.1191300000000_dp,
     >  -0.2062660000000_dp,   0.6673880000000_dp,   0.2718070000000_dp,
     >   0.4405230000000_dp,   0.1179360000000_dp,   0.7250570000000_dp,
     >  -0.4002600000000_dp,   0.5441920000000_dp,  -0.0666994000000_dp,
     >   0.3716500000000_dp
     >      /), (/20,20/))

      real(kind=dp), parameter, dimension(20) ::
     >  b4 = (/ 
     >  -0.0446875000000_dp,  -0.4246940000000_dp,  -0.2270150000000_dp,
     >   0.3740400000000_dp,   0.7332760000000_dp,  -0.1108520000000_dp,
     >  -0.7114290000000_dp,   1.2662100000000_dp,   0.8961380000000_dp,
     >   0.2500370000000_dp,   0.2555510000000_dp,   0.2565220000000_dp,
     >  -0.0361706000000_dp,   0.0528172000000_dp,   0.3192980000000_dp,
     >  -0.2791350000000_dp,   0.2675700000000_dp,   0.4969490000000_dp,
     >  -0.7693630000000_dp,   0.8420010000000_dp
     >      /)

      real(kind=dp), parameter, dimension(20,20) ::
     >  W4 = reshape( (/ 
     >   0.3681770000000_dp,  -1.3654800000000_dp,   1.0482600000000_dp,
     >   0.0402462000000_dp,  -0.0809720000000_dp,  -0.7810610000000_dp,
     >  -1.4281300000000_dp,  -0.0913602000000_dp,   0.0461368000000_dp,
     >  -0.1491430000000_dp,   0.0995101000000_dp,   0.6842060000000_dp,
     >   0.3704420000000_dp,   1.0077300000000_dp,  -0.0513708000000_dp,
     >  -1.8609700000000_dp,   0.5215760000000_dp,  -0.9106060000000_dp,
     >  -0.3543550000000_dp,   0.8150030000000_dp,   0.1848290000000_dp,
     >   0.5765930000000_dp,  -1.1225900000000_dp,   0.4251290000000_dp,
     >  -0.2245820000000_dp,   1.2114400000000_dp,   1.3085200000000_dp,
     >  -0.2800920000000_dp,  -0.2611790000000_dp,   0.1196320000000_dp,
     >  -0.1500170000000_dp,  -1.0534100000000_dp,  -0.6477350000000_dp,
     >  -0.5347540000000_dp,   0.5472240000000_dp,   0.5352840000000_dp,
     >  -0.7499990000000_dp,   1.1318900000000_dp,  -0.1941270000000_dp,
     >  -0.4370420000000_dp,  -0.5602270000000_dp,   2.6101100000000_dp,
     >  -4.5945500000000_dp,   0.9468040000000_dp,  -0.0751428000000_dp,
     >   0.8100330000000_dp,   2.3571900000000_dp,  -0.0143629000000_dp,
     >   0.4134260000000_dp,   0.4054230000000_dp,   1.1348800000000_dp,
     >  -0.3768360000000_dp,  -0.5654470000000_dp,  -1.3001600000000_dp,
     >   0.2701600000000_dp,   1.8639600000000_dp,   0.0745336000000_dp,
     >   3.5620000000000_dp,   0.1327120000000_dp,  -3.4882700000000_dp,
     >  -0.1151850000000_dp,   0.0431694000000_dp,   0.5289690000000_dp,
     >   0.0814686000000_dp,   0.4992320000000_dp,  -0.2391230000000_dp,
     >  -0.3986030000000_dp,   0.7742660000000_dp,  -0.1515230000000_dp,
     >  -0.0644823000000_dp,   0.3392280000000_dp,   0.3630100000000_dp,
     >  -0.2065170000000_dp,   0.0768753000000_dp,  -0.2323180000000_dp,
     >  -0.0981481000000_dp,   0.1549670000000_dp,   0.7209150000000_dp,
     >   0.0681167000000_dp,   0.3372710000000_dp,  -0.4620340000000_dp,
     >   1.0318000000000_dp,  -2.1872000000000_dp,   0.0777344000000_dp,
     >   0.3206220000000_dp,   0.7695740000000_dp,   1.2825200000000_dp,
     >   0.4043690000000_dp,  -0.2868660000000_dp,   0.5082750000000_dp,
     >   0.3289900000000_dp,  -0.3483470000000_dp,  -0.6048420000000_dp,
     >  -1.0739700000000_dp,   0.4129650000000_dp,   1.3090800000000_dp,
     >  -0.3393940000000_dp,   1.6599000000000_dp,   0.3756140000000_dp,
     >  -1.4716600000000_dp,   0.1645140000000_dp,  -1.6317200000000_dp,
     >   2.1890000000000_dp,   0.1351540000000_dp,  -0.0123093000000_dp,
     >  -0.8386050000000_dp,  -2.3959700000000_dp,  -0.1232090000000_dp,
     >   0.0201091000000_dp,   0.0121008000000_dp,   0.1138070000000_dp,
     >   0.1611150000000_dp,   0.3775430000000_dp,   0.7651660000000_dp,
     >  -0.0579404000000_dp,  -1.9279500000000_dp,   0.3282660000000_dp,
     >  -1.6071700000000_dp,   0.0074790400000_dp,   1.7660800000000_dp,
     >  -0.1202860000000_dp,   0.7656990000000_dp,  -0.1109630000000_dp,
     >  -0.1187500000000_dp,  -0.9001140000000_dp,   0.2916690000000_dp,
     >   0.4095050000000_dp,   0.1080590000000_dp,  -0.8973460000000_dp,
     >   0.5325320000000_dp,  -0.5027450000000_dp,  -0.9979390000000_dp,
     >  -0.5023390000000_dp,  -0.1041990000000_dp,  -0.3766520000000_dp,
     >   0.4980770000000_dp,  -0.7329150000000_dp,   0.3467190000000_dp,
     >   1.0142400000000_dp,   0.0914130000000_dp,   0.5368350000000_dp,
     >   0.2101030000000_dp,  -0.2922420000000_dp,   0.5342860000000_dp,
     >  -0.1611220000000_dp,  -0.5151930000000_dp,  -0.2180050000000_dp,
     >  -0.0572768000000_dp,   0.4097450000000_dp,   0.2952060000000_dp,
     >   0.7333660000000_dp,   0.5522300000000_dp,   0.4689820000000_dp,
     >   0.5701070000000_dp,   0.1789790000000_dp,  -0.6043320000000_dp,
     >   0.5951000000000_dp,   0.0066079300000_dp,  -0.2695290000000_dp,
     >   0.1196160000000_dp,  -0.1477270000000_dp,   0.0680111000000_dp,
     >  -1.0474000000000_dp,   0.5931870000000_dp,  -0.1699430000000_dp,
     >   0.7123690000000_dp,   0.3930110000000_dp,  -0.0634108000000_dp,
     >   0.2700300000000_dp,  -0.5808850000000_dp,  -0.4442510000000_dp,
     >  -0.5466270000000_dp,  -0.6984310000000_dp,  -0.6033740000000_dp,
     >  -0.8169200000000_dp,   0.8858340000000_dp,  -0.2699940000000_dp,
     >   0.8134250000000_dp,   0.5020720000000_dp,  -0.3153560000000_dp,
     >  -0.3760150000000_dp,   1.7588900000000_dp,  -2.9859800000000_dp,
     >   0.2587560000000_dp,  -0.3196010000000_dp,   0.7349340000000_dp,
     >   1.9273100000000_dp,  -0.6661700000000_dp,   0.1067060000000_dp,
     >   0.1569330000000_dp,   0.2407840000000_dp,  -0.5160840000000_dp,
     >  -0.3304990000000_dp,  -1.1685800000000_dp,  -0.2925350000000_dp,
     >   1.3840600000000_dp,  -0.0026037600000_dp,   1.9700500000000_dp,
     >   0.2042830000000_dp,  -2.5952100000000_dp,  -1.2646700000000_dp,
     >   2.0800500000000_dp,  -2.0485800000000_dp,   0.1282520000000_dp,
     >   0.0185757000000_dp,   0.6913800000000_dp,   1.8641500000000_dp,
     >  -0.0063702500000_dp,   0.2250120000000_dp,   0.7640860000000_dp,
     >   0.4859410000000_dp,  -0.1719140000000_dp,  -0.6742390000000_dp,
     >  -1.1922500000000_dp,  -0.3178010000000_dp,   2.4244500000000_dp,
     >  -0.4213350000000_dp,   2.4931800000000_dp,  -0.1181320000000_dp,
     >  -1.8007500000000_dp,  -0.0296414000000_dp,  -1.0619800000000_dp,
     >  -0.7986890000000_dp,   0.2138200000000_dp,   0.6207440000000_dp,
     >   0.9995780000000_dp,   1.0449800000000_dp,  -0.2181430000000_dp,
     >   0.2859060000000_dp,  -0.5087760000000_dp,   0.3442760000000_dp,
     >  -0.2535400000000_dp,  -0.6297460000000_dp,  -1.3825800000000_dp,
     >   0.2803470000000_dp,   0.7760680000000_dp,  -0.0028701900000_dp,
     >   0.8039470000000_dp,  -0.6302610000000_dp,  -0.5345020000000_dp,
     >   0.7720140000000_dp,  -0.7637340000000_dp,   0.2992290000000_dp,
     >   0.2135340000000_dp,  -0.1097810000000_dp,  -0.2237550000000_dp,
     >  -1.0639700000000_dp,  -0.4985150000000_dp,   0.2043030000000_dp,
     >   0.2178020000000_dp,   0.0354621000000_dp,  -0.0174415000000_dp,
     >  -0.0289828000000_dp,   0.7444310000000_dp,  -1.4205500000000_dp,
     >  -1.3503500000000_dp,  -0.1910030000000_dp,  -0.7492280000000_dp,
     >   0.0493873000000_dp,   0.4408970000000_dp,  -0.0545817000000_dp,
     >  -0.5191960000000_dp,   0.3576610000000_dp,  -0.2775740000000_dp,
     >  -0.1711430000000_dp,  -0.6769510000000_dp,  -0.9119910000000_dp,
     >   0.1654590000000_dp,  -0.1729500000000_dp,  -0.4177040000000_dp,
     >  -0.1566370000000_dp,   0.7246580000000_dp,   0.8199890000000_dp,
     >   0.4626800000000_dp,   1.0415100000000_dp,  -1.0320900000000_dp,
     >   0.5437900000000_dp,  -0.4919590000000_dp,  -0.3933160000000_dp,
     >   0.2248830000000_dp,  -0.3257230000000_dp,  -0.1616270000000_dp,
     >  -1.6425300000000_dp,  -0.1342060000000_dp,   0.0898834000000_dp,
     >   0.6904480000000_dp,   0.7570480000000_dp,  -0.2081310000000_dp,
     >  -0.2169510000000_dp,   0.3366160000000_dp,  -0.0018988400000_dp,
     >  -0.3405000000000_dp,  -0.4400760000000_dp,  -0.9433380000000_dp,
     >   0.1553210000000_dp,   0.4027180000000_dp,  -0.2780020000000_dp,
     >   1.0342800000000_dp,   0.0495943000000_dp,  -0.9890740000000_dp,
     >   0.5093500000000_dp,  -1.8133700000000_dp,   3.0327100000000_dp,
     >  -0.3415200000000_dp,   0.3382860000000_dp,  -0.7673450000000_dp,
     >  -2.5945400000000_dp,  -0.4237980000000_dp,   0.2201670000000_dp,
     >  -0.6320410000000_dp,  -0.6463880000000_dp,   0.0376426000000_dp,
     >   0.2686440000000_dp,   1.2744800000000_dp,   0.1967930000000_dp,
     >  -1.7053700000000_dp,   0.2516490000000_dp,  -1.8941600000000_dp,
     >  -0.4714910000000_dp,   2.6667000000000_dp,   0.2449910000000_dp,
     >   1.0730900000000_dp,  -2.5076500000000_dp,   0.4438000000000_dp,
     >  -0.0195254000000_dp,   0.6941290000000_dp,   1.7134300000000_dp,
     >   0.0144911000000_dp,   0.0588716000000_dp,   0.3822930000000_dp,
     >   0.3622740000000_dp,  -0.2965120000000_dp,  -0.6314390000000_dp,
     >  -1.0890900000000_dp,  -0.3780700000000_dp,   1.4680400000000_dp,
     >  -0.2754070000000_dp,   1.7396100000000_dp,   0.2977360000000_dp,
     >  -1.3483400000000_dp,   0.9991560000000_dp,  -0.7995070000000_dp,
     >   1.7119700000000_dp,   0.4200380000000_dp,  -0.0111848000000_dp,
     >  -0.6826420000000_dp,  -1.3363300000000_dp,   0.3857000000000_dp,
     >   0.2602610000000_dp,   0.3349950000000_dp,  -0.1825800000000_dp,
     >   0.2751210000000_dp,   0.3147790000000_dp,   1.4837300000000_dp,
     >   0.2518050000000_dp,  -1.2560400000000_dp,   0.1413040000000_dp,
     >  -1.9441000000000_dp,   0.8146370000000_dp,   0.9113900000000_dp,
     >   0.2724330000000_dp,  -1.5400500000000_dp,   1.1407300000000_dp,
     >  -0.0809263000000_dp,   0.2763470000000_dp,  -0.3744970000000_dp,
     >  -1.1657600000000_dp,   0.1044560000000_dp,   0.1533430000000_dp,
     >   0.0890387000000_dp,   0.2817500000000_dp,   1.0170700000000_dp,
     >  -0.0980648000000_dp,   0.4439660000000_dp,   0.2473160000000_dp,
     >  -0.9038530000000_dp,   0.2558480000000_dp,  -0.3793440000000_dp,
     >  -0.2361900000000_dp,   0.9770690000000_dp,   0.3118710000000_dp,
     >  -2.0215800000000_dp,   1.4802600000000_dp,  -0.3302370000000_dp,
     >   0.5870760000000_dp,  -0.7347180000000_dp,  -2.1770500000000_dp,
     >   0.5654870000000_dp,   0.7470520000000_dp,  -0.1667580000000_dp,
     >  -0.5796720000000_dp,   0.2749070000000_dp,   0.7220440000000_dp,
     >   1.3101000000000_dp,   0.2946910000000_dp,  -2.6166500000000_dp,
     >   0.1736810000000_dp,  -1.5748100000000_dp,  -0.0360835000000_dp,
     >   1.5336000000000_dp
     >      /), (/20,20/))

      real(kind=dp), parameter ::
     >  b5 =   -0.4182790000000_dp

      real(kind=dp), parameter, dimension(20) ::
     >  W5 = (/ 
     >   0.3377520000000_dp,   0.3740840000000_dp,  -0.5932940000000_dp,
     >  -0.3614160000000_dp,  -0.1552180000000_dp,   0.4646470000000_dp,
     >   0.6654670000000_dp,  -0.2658030000000_dp,  -0.2026280000000_dp,
     >  -0.2158520000000_dp,   0.2409480000000_dp,  -0.2813920000000_dp,
     >  -0.3152220000000_dp,  -0.5652230000000_dp,  -0.1000420000000_dp,
     >  -0.6429410000000_dp,  -0.2900930000000_dp,  -0.5437890000000_dp,
     >   0.1843680000000_dp,   0.5288490000000_dp
     >      /)

C======================================================================|