%COPPER2  Black-copper-white colormap
%
% Examples:
%   map = copper2;
%   map = copper2(len);
%   B = copper2(A);
%   B = copper2(A, lims);
%
% A brighter variation of MATLAB's copper colormap, which also continues on
% to white. This colormap has been designed to achieve perceptual
% uniformity.
%
% The function can additionally be used to convert a real-valued array into
% a truecolor array using the colormap.
%
% IN:
%   len - Scalar length of the output colormap. If len == Inf the concise
%         table is returned. Default: len = size(get(gcf, 'Colormap'), 1);
%   A - Non-scalar numeric array of real values to be converted into
%       truecolor.
%   lims - 1x2 array of saturation limits to be used on A. Default:
%          [min(A(:)) max(A(:))].
%
% OUT:
%   map - (len)x3 colormap table.
%   B - size(A)x3 truecolor array.

% Copyright: Oliver Woodford, 2009, 2015

function map = copper2(varargin)
map = [     0         0         0
       0.0024    0.0069    0.0018
       0.0049    0.0139    0.0035
       0.0074    0.0208    0.0053
       0.0098    0.0277    0.0070
       0.0126    0.0346    0.0087
       0.0151    0.0416    0.0103
       0.0178    0.0478    0.0122
       0.0205    0.0537    0.0138
       0.0233    0.0592    0.0153
       0.0262    0.0642    0.0167
       0.0291    0.0690    0.0181
       0.0322    0.0735    0.0195
       0.0357    0.0777    0.0207
       0.0390    0.0819    0.0219
       0.0424    0.0857    0.0230
       0.0457    0.0895    0.0241
       0.0492    0.0930    0.0250
       0.0525    0.0964    0.0258
       0.0560    0.0998    0.0265
       0.0594    0.1030    0.0272
       0.0628    0.1061    0.0278
       0.0659    0.1093    0.0282
       0.0691    0.1125    0.0285
       0.0719    0.1157    0.0287
       0.0747    0.1190    0.0288
       0.0773    0.1223    0.0287
       0.0799    0.1257    0.0285
       0.0827    0.1290    0.0282
       0.0855    0.1324    0.0279
       0.0885    0.1357    0.0274
       0.0917    0.1390    0.0268
       0.0951    0.1423    0.0262
       0.0986    0.1456    0.0256
       0.1023    0.1488    0.0248
       0.1062    0.1520    0.0240
       0.1103    0.1552    0.0233
       0.1144    0.1583    0.0224
       0.1189    0.1615    0.0216
       0.1234    0.1646    0.0208
       0.1281    0.1676    0.0200
       0.1329    0.1707    0.0192
       0.1379    0.1737    0.0185
       0.1431    0.1766    0.0179
       0.1483    0.1796    0.0173
       0.1537    0.1825    0.0169
       0.1593    0.1854    0.0165
       0.1649    0.1882    0.0162
       0.1707    0.1910    0.0161
       0.1766    0.1937    0.0160
       0.1826    0.1965    0.0162
       0.1887    0.1991    0.0165
       0.1949    0.2018    0.0170
       0.2012    0.2044    0.0177
       0.2076    0.2070    0.0185
       0.2141    0.2095    0.0195
       0.2206    0.2121    0.0207
       0.2272    0.2145    0.0220
       0.2338    0.2170    0.0236
       0.2405    0.2194    0.0255
       0.2473    0.2218    0.0275
       0.2541    0.2242    0.0297
       0.2609    0.2265    0.0321
       0.2678    0.2289    0.0347
       0.2747    0.2312    0.0376
       0.2817    0.2334    0.0405
       0.2886    0.2357    0.0434
       0.2956    0.2379    0.0468
       0.3026    0.2401    0.0502
       0.3096    0.2423    0.0538
       0.3165    0.2445    0.0582
       0.3235    0.2467    0.0617
       0.3305    0.2489    0.0657
       0.3375    0.2511    0.0698
       0.3444    0.2532    0.0740
       0.3514    0.2554    0.0782
       0.3583    0.2576    0.0827
       0.3652    0.2597    0.0872
       0.3721    0.2619    0.0919
       0.3790    0.2640    0.0965
       0.3858    0.2662    0.1013
       0.3927    0.2683    0.1062
       0.3995    0.2704    0.1110
       0.4064    0.2726    0.1159
       0.4132    0.2747    0.1208
       0.4200    0.2768    0.1259
       0.4268    0.2790    0.1309
       0.4336    0.2811    0.1358
       0.4405    0.2832    0.1406
       0.4473    0.2853    0.1456
       0.4542    0.2874    0.1507
       0.4610    0.2895    0.1557
       0.4678    0.2916    0.1608
       0.4746    0.2937    0.1658
       0.4815    0.2958    0.1708
       0.4884    0.2979    0.1758
       0.4953    0.2999    0.1808
       0.5022    0.3019    0.1856
       0.5092    0.3040    0.1904
       0.5162    0.3060    0.1952
       0.5232    0.3079    0.1999
       0.5303    0.3099    0.2047
       0.5373    0.3119    0.2095
       0.5445    0.3138    0.2140
       0.5515    0.3158    0.2188
       0.5586    0.3177    0.2235
       0.5658    0.3197    0.2282
       0.5730    0.3215    0.2327
       0.5802    0.3237    0.2342
       0.5873    0.3258    0.2356
       0.5945    0.3279    0.2370
       0.6017    0.3300    0.2384
       0.6089    0.3321    0.2397
       0.6161    0.3343    0.2410
       0.6232    0.3364    0.2422
       0.6304    0.3386    0.2434
       0.6377    0.3406    0.2452
       0.6451    0.3426    0.2471
       0.6525    0.3446    0.2489
       0.6599    0.3466    0.2506
       0.6672    0.3486    0.2523
       0.6746    0.3506    0.2539
       0.6819    0.3527    0.2554
       0.6893    0.3547    0.2569
       0.6956    0.3576    0.2568
       0.7014    0.3608    0.2563
       0.7072    0.3641    0.2557
       0.7129    0.3674    0.2551
       0.7185    0.3708    0.2544
       0.7240    0.3743    0.2537
       0.7295    0.3779    0.2530
       0.7349    0.3815    0.2522
       0.7404    0.3850    0.2515
       0.7461    0.3884    0.2509
       0.7518    0.3918    0.2502
       0.7574    0.3953    0.2495
       0.7628    0.3989    0.2487
       0.7682    0.4025    0.2480
       0.7736    0.4062    0.2471
       0.7788    0.4100    0.2463
       0.7837    0.4140    0.2454
       0.7873    0.4189    0.2449
       0.7907    0.4238    0.2445
       0.7941    0.4288    0.2441
       0.7974    0.4338    0.2438
       0.8006    0.4389    0.2436
       0.8037    0.4441    0.2434
       0.8067    0.4492    0.2433
       0.8097    0.4544    0.2433
       0.8126    0.4597    0.2433
       0.8154    0.4649    0.2435
       0.8181    0.4702    0.2436
       0.8208    0.4755    0.2439
       0.8234    0.4809    0.2443
       0.8259    0.4862    0.2447
       0.8284    0.4917    0.2452
       0.8308    0.4971    0.2458
       0.8324    0.5029    0.2480
       0.8338    0.5088    0.2509
       0.8352    0.5146    0.2539
       0.8365    0.5205    0.2571
       0.8378    0.5263    0.2605
       0.8390    0.5322    0.2640
       0.8403    0.5380    0.2677
       0.8414    0.5438    0.2717
       0.8426    0.5497    0.2759
       0.8437    0.5555    0.2804
       0.8447    0.5612    0.2851
       0.8458    0.5670    0.2900
       0.8469    0.5728    0.2950
       0.8479    0.5785    0.3003
       0.8490    0.5842    0.3057
       0.8501    0.5898    0.3114
       0.8512    0.5955    0.3172
       0.8523    0.6011    0.3231
       0.8534    0.6066    0.3293
       0.8546    0.6122    0.3356
       0.8558    0.6176    0.3421
       0.8571    0.6231    0.3487
       0.8584    0.6285    0.3556
       0.8598    0.6339    0.3625
       0.8612    0.6392    0.3697
       0.8627    0.6444    0.3770
       0.8643    0.6496    0.3845
       0.8660    0.6548    0.3921
       0.8679    0.6598    0.3998
       0.8697    0.6649    0.4076
       0.8718    0.6699    0.4156
       0.8739    0.6748    0.4236
       0.8761    0.6796    0.4317
       0.8784    0.6845    0.4399
       0.8807    0.6893    0.4481
       0.8832    0.6940    0.4564
       0.8856    0.6988    0.4647
       0.8881    0.7035    0.4730
       0.8906    0.7082    0.4813
       0.8932    0.7129    0.4897
       0.8957    0.7176    0.4980
       0.8982    0.7223    0.5064
       0.9006    0.7270    0.5148
       0.9031    0.7318    0.5231
       0.9055    0.7365    0.5315
       0.9080    0.7412    0.5399
       0.9104    0.7459    0.5483
       0.9127    0.7507    0.5567
       0.9151    0.7554    0.5651
       0.9174    0.7602    0.5735
       0.9197    0.7649    0.5819
       0.9220    0.7697    0.5904
       0.9243    0.7744    0.5988
       0.9265    0.7792    0.6073
       0.9288    0.7840    0.6157
       0.9310    0.7888    0.6242
       0.9331    0.7936    0.6327
       0.9353    0.7984    0.6411
       0.9374    0.8032    0.6497
       0.9395    0.8080    0.6582
       0.9415    0.8128    0.6667
       0.9436    0.8176    0.6752
       0.9456    0.8224    0.6838
       0.9476    0.8273    0.6923
       0.9496    0.8321    0.7009
       0.9515    0.8370    0.7095
       0.9534    0.8418    0.7181
       0.9553    0.8467    0.7267
       0.9572    0.8515    0.7353
       0.9590    0.8564    0.7439
       0.9608    0.8613    0.7526
       0.9626    0.8661    0.7612
       0.9643    0.8710    0.7699
       0.9661    0.8759    0.7786
       0.9678    0.8808    0.7873
       0.9695    0.8857    0.7960
       0.9711    0.8906    0.8047
       0.9727    0.8955    0.8134
       0.9743    0.9005    0.8222
       0.9759    0.9054    0.8310
       0.9774    0.9103    0.8397
       0.9789    0.9152    0.8485
       0.9804    0.9202    0.8573
       0.9819    0.9251    0.8661
       0.9833    0.9301    0.8750
       0.9847    0.9350    0.8838
       0.9860    0.9400    0.8927
       0.9874    0.9450    0.9015
       0.9886    0.9500    0.9104
       0.9899    0.9549    0.9193
       0.9912    0.9599    0.9282
       0.9924    0.9649    0.9372
       0.9935    0.9699    0.9461
       0.9947    0.9749    0.9551
       0.9958    0.9799    0.9640
       0.9969    0.9849    0.9730
       0.9980    0.9900    0.9820
       0.9990    0.9950    0.9910
            1         1         1];
map = colormap_helper(map, varargin{:});