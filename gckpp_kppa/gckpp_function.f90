!----------------------- BEGIN gckpp_function.f90 BEGIN ----------------------
! @file gckpp_function.f90                                                    
! @author yifanwu                                                             
! @date 2014-04-08 15:56:45.403181                                            
! @brief The ODE function of the chemical model                               
!                                                                             
! The ODE function of the chemical model                                      
!                                                                             
! This file was generated by Kppa: http://www.paratools.com/Kppa              
!-----------------------------------------------------------------------------


MODULE gckpp_function

  USE gckpp_parameters
  USE gckpp_sparse

  IMPLICIT NONE



  CONTAINS

!------------------------------------ Fun ------------------------------------
! The ODE function of the chemical model                                      
!                                                                             
! @param[in]     var    Variable species concentrations                       
! @param[in]     fix    Fixed species concentrations                          
! @param[in]     rct    Reaction rates                                        
! @param[out]    vardot The ODE function                                      
!-----------------------------------------------------------------------------
  SUBROUTINE Fun(var, fix, rct, vardot)
    IMPLICIT NONE

    REAL(8), INTENT(IN) :: var(NVAR)
    REAL(8), INTENT(IN) :: fix(NFIX)
    REAL(8), INTENT(IN) :: rct(NREACT)
    REAL(8), INTENT(OUT) :: vardot(NVAR)

    REAL(8) :: r0
    REAL(8) :: r1
    REAL(8) :: r2
    REAL(8) :: r3
    REAL(8) :: r4
    REAL(8) :: r5
    REAL(8) :: r6
    REAL(8) :: r7
    REAL(8) :: r8
    REAL(8) :: r9
    REAL(8) :: r10
    REAL(8) :: r11
    REAL(8) :: r12
    REAL(8) :: r13
    REAL(8) :: r14
    REAL(8) :: r15
    REAL(8) :: r16
    REAL(8) :: r17
    REAL(8) :: r18
    REAL(8) :: r19
    REAL(8) :: r20
    REAL(8) :: r21
    REAL(8) :: r22
    REAL(8) :: r23
    REAL(8) :: r24
    REAL(8) :: r25
    REAL(8) :: r26
    REAL(8) :: r27
    REAL(8) :: r28
    REAL(8) :: r29
    REAL(8) :: r30
    REAL(8) :: r31
    REAL(8) :: r32
    REAL(8) :: r33
    REAL(8) :: r34
    REAL(8) :: r35
    REAL(8) :: r36
    REAL(8) :: r37
    REAL(8) :: r38
    REAL(8) :: r39
    REAL(8) :: r40
    REAL(8) :: r41
    REAL(8) :: r42
    REAL(8) :: r43
    REAL(8) :: r44
    REAL(8) :: r45
    REAL(8) :: r46
    REAL(8) :: r47
    REAL(8) :: r48
    REAL(8) :: r49
    REAL(8) :: r50
    REAL(8) :: r51
    REAL(8) :: r52
    REAL(8) :: r53
    REAL(8) :: r54
    REAL(8) :: r55
    REAL(8) :: r56
    REAL(8) :: r57
    REAL(8) :: r58
    REAL(8) :: r59
    REAL(8) :: r60
    REAL(8) :: r61
    REAL(8) :: r62
    REAL(8) :: r63
    REAL(8) :: r64
    REAL(8) :: r65
    REAL(8) :: r66
    REAL(8) :: r67
    REAL(8) :: r68
    REAL(8) :: r69
    REAL(8) :: r70
    REAL(8) :: r71
    REAL(8) :: r72
    REAL(8) :: r73
    REAL(8) :: r74
    REAL(8) :: r75
    REAL(8) :: r76
    REAL(8) :: r77
    REAL(8) :: r78
    REAL(8) :: r79
    REAL(8) :: r80
    REAL(8) :: r81
    REAL(8) :: r82
    REAL(8) :: r83
    REAL(8) :: r84
    REAL(8) :: r85
    REAL(8) :: r86
    REAL(8) :: r87
    REAL(8) :: r88
    REAL(8) :: r89
    REAL(8) :: r90
    REAL(8) :: r91
    REAL(8) :: r92
    REAL(8) :: r93
    REAL(8) :: r94
    REAL(8) :: r95
    REAL(8) :: r96
    REAL(8) :: r97
    REAL(8) :: r98
    REAL(8) :: r99
    REAL(8) :: r100
    REAL(8) :: r101
    REAL(8) :: r102
    REAL(8) :: r103
    REAL(8) :: r104
    REAL(8) :: r105
    REAL(8) :: r106
    REAL(8) :: r107
    REAL(8) :: r108
    REAL(8) :: r109
    REAL(8) :: r110
    REAL(8) :: r111
    REAL(8) :: r112
    REAL(8) :: r113
    REAL(8) :: r114
    REAL(8) :: r115
    REAL(8) :: r116
    REAL(8) :: r117
    REAL(8) :: r118
    REAL(8) :: r119
    REAL(8) :: r120
    REAL(8) :: r121
    REAL(8) :: r122
    REAL(8) :: r123
    REAL(8) :: r124
    REAL(8) :: r125
    REAL(8) :: r126
    REAL(8) :: r127
    REAL(8) :: r128
    REAL(8) :: r129
    REAL(8) :: r130
    REAL(8) :: r131
    REAL(8) :: r132
    REAL(8) :: r133
    REAL(8) :: r134
    REAL(8) :: r135
    REAL(8) :: r136
    REAL(8) :: r137
    REAL(8) :: r138
    REAL(8) :: r139
    REAL(8) :: r140
    REAL(8) :: r141
    REAL(8) :: r142
    REAL(8) :: r143
    REAL(8) :: r144
    REAL(8) :: r145
    REAL(8) :: r146
    REAL(8) :: r147
    REAL(8) :: r148
    REAL(8) :: r149
    REAL(8) :: r150
    REAL(8) :: r151
    REAL(8) :: r152
    REAL(8) :: r153
    REAL(8) :: r154
    REAL(8) :: r155
    REAL(8) :: r156
    REAL(8) :: r157
    REAL(8) :: r158
    REAL(8) :: r159
    REAL(8) :: r160
    REAL(8) :: r161
    REAL(8) :: r162
    REAL(8) :: r163
    REAL(8) :: r164
    REAL(8) :: r165
    REAL(8) :: r166
    REAL(8) :: r167
    REAL(8) :: r168
    REAL(8) :: r169
    REAL(8) :: r170
    REAL(8) :: r171
    REAL(8) :: r172
    REAL(8) :: r173
    REAL(8) :: r174
    REAL(8) :: r175
    REAL(8) :: r176
    REAL(8) :: r177
    REAL(8) :: r178
    REAL(8) :: r179
    REAL(8) :: r180
    REAL(8) :: r181
    REAL(8) :: r182
    REAL(8) :: r183
    REAL(8) :: r184
    REAL(8) :: r185
    REAL(8) :: r186
    REAL(8) :: r187
    REAL(8) :: r188
    REAL(8) :: r189
    REAL(8) :: r190
    REAL(8) :: r191
    REAL(8) :: r192
    REAL(8) :: r193
    REAL(8) :: r194
    REAL(8) :: r195
    REAL(8) :: r196
    REAL(8) :: r197
    REAL(8) :: r198
    REAL(8) :: r199
    REAL(8) :: r200
    REAL(8) :: r201
    REAL(8) :: r202
    REAL(8) :: r203
    REAL(8) :: r204
    REAL(8) :: r205
    REAL(8) :: r206
    REAL(8) :: r207
    REAL(8) :: r208
    REAL(8) :: r209
    REAL(8) :: r210
    REAL(8) :: r211
    REAL(8) :: r212
    REAL(8) :: r213
    REAL(8) :: r214
    REAL(8) :: r215
    REAL(8) :: r216
    REAL(8) :: r217
    REAL(8) :: r218
    REAL(8) :: r219
    REAL(8) :: r220
    REAL(8) :: r221
    REAL(8) :: r222
    REAL(8) :: r223
    REAL(8) :: r224
    REAL(8) :: r225
    REAL(8) :: r226
    REAL(8) :: r227
    REAL(8) :: r228
    REAL(8) :: r229
    REAL(8) :: r230
    REAL(8) :: r231
    REAL(8) :: r232
    REAL(8) :: r233
    REAL(8) :: r234
    REAL(8) :: r235
    REAL(8) :: r236
    REAL(8) :: r237
    REAL(8) :: r238
    REAL(8) :: r239
    REAL(8) :: r240
    REAL(8) :: r241
    REAL(8) :: r242
    REAL(8) :: r243
    REAL(8) :: r244
    REAL(8) :: r245
    REAL(8) :: r246
    REAL(8) :: r247
    REAL(8) :: r248
    REAL(8) :: r249
    REAL(8) :: r250
    REAL(8) :: r251
    REAL(8) :: r252
    REAL(8) :: r253
    REAL(8) :: r254
    REAL(8) :: r255
    REAL(8) :: r256
    REAL(8) :: r257
    REAL(8) :: r258
    REAL(8) :: r259
    REAL(8) :: r260
    REAL(8) :: r261
    REAL(8) :: r262
    REAL(8) :: r263
    REAL(8) :: r264
    REAL(8) :: r265
    REAL(8) :: r266
    REAL(8) :: r267
    REAL(8) :: r268
    REAL(8) :: r269
    REAL(8) :: r270
    REAL(8) :: r271
    REAL(8) :: r272
    REAL(8) :: r273
    REAL(8) :: r274
    REAL(8) :: r275
    REAL(8) :: r276
    REAL(8) :: r277
    REAL(8) :: r278
    REAL(8) :: r279
    REAL(8) :: r280
    REAL(8) :: r281
    REAL(8) :: r282
    REAL(8) :: r283
    REAL(8) :: r284
    REAL(8) :: r285
    REAL(8) :: r286
    REAL(8) :: r287
    REAL(8) :: r288
    REAL(8) :: r289
    REAL(8) :: r290
    REAL(8) :: r291
    REAL(8) :: r292
    REAL(8) :: r293
    REAL(8) :: r294
    REAL(8) :: r295
    REAL(8) :: r296
    REAL(8) :: r297
    REAL(8) :: r298
    REAL(8) :: r299
    REAL(8) :: r300
    REAL(8) :: r301
    REAL(8) :: r302
    REAL(8) :: r303
    REAL(8) :: r304
    REAL(8) :: r305
    REAL(8) :: r306
    REAL(8) :: r307
    REAL(8) :: r308
    REAL(8) :: r309
    REAL(8) :: r310
    REAL(8) :: r311
    REAL(8) :: r312
    REAL(8) :: r313
    REAL(8) :: r314
    REAL(8) :: r315
    REAL(8) :: r316
    REAL(8) :: r317
    REAL(8) :: r318
    REAL(8) :: r319
    REAL(8) :: r320
    REAL(8) :: r321
    REAL(8) :: r322
    REAL(8) :: r323
    REAL(8) :: r324
    REAL(8) :: r325
    REAL(8) :: r326
    REAL(8) :: r327
    REAL(8) :: r328
    REAL(8) :: r329
    REAL(8) :: r330

    r0 = rct(1)*var(83)*var(90)
    r1 = rct(2)*var(89)*var(90)
    r2 = rct(3)*var(84)*var(90)
    r3 = rct(4)*var(87)*var(90)
    r4 = rct(5)*var(90)*var(92)
    r5 = rct(6)*var(89)**2
    r6 = rct(7)*var(89)**2
    r7 = rct(8)*var(84)*var(89)
    r8 = rct(9)*var(11)*var(89)
    r9 = rct(10)*var(83)*var(84)
    r10 = rct(11)*var(84)**2
    r11 = fix(6)*rct(12)*var(89)
    r12 = rct(13)*var(48)*var(89)
    r13 = fix(2)*rct(14)*var(89)
    r14 = rct(15)*var(83)*var(92)
    r15 = rct(16)*var(84)*var(92)
    r16 = rct(17)*var(92)**2
    r17 = rct(18)*var(92)**2
    r18 = rct(19)*var(22)*var(89)
    r19 = rct(20)*var(22)*var(89)
    r20 = rct(21)*var(72)*var(89)
    r21 = rct(22)*var(87)*var(89)
    r22 = rct(23)*var(54)*var(89)
    r23 = rct(24)*var(83)*var(89)
    r24 = rct(25)*var(17)*var(89)
    r25 = rct(26)*var(84)*var(87)
    r26 = rct(27)*var(19)
    r27 = rct(28)*var(19)*var(89)
    r28 = rct(29)*var(84)*var(85)
    r29 = rct(30)*var(83)*var(85)
    r30 = rct(31)*var(85)*var(89)
    r31 = rct(32)*var(85)*var(87)
    r32 = rct(33)*var(18)
    r33 = fix(9)*rct(34)*var(89)
    r34 = rct(35)*var(85)*var(87)
    r35 = rct(36)*var(72)*var(85)
    r36 = rct(37)*var(80)*var(89)
    r37 = rct(38)*var(80)*var(85)
    r38 = rct(39)*var(87)*var(88)
    r39 = rct(40)*var(15)
    r40 = rct(41)*var(83)*var(88)
    r41 = rct(42)*var(40)*var(89)
    r42 = rct(43)*var(82)*var(83)
    r43 = rct(44)*var(30)*var(89)
    r44 = rct(45)*var(30)*var(89)
    r45 = rct(46)*var(55)*var(83)
    r46 = rct(47)*var(65)*var(83)
    r47 = rct(48)*var(16)*var(89)
    r48 = rct(49)*var(74)*var(83)
    r49 = rct(50)*var(74)*var(83)
    r50 = rct(51)*var(56)*var(83)
    r51 = rct(52)*var(67)*var(83)
    r52 = rct(53)*var(70)*var(83)
    r53 = rct(54)*var(78)*var(83)
    r54 = rct(55)*var(58)*var(83)
    r55 = rct(56)*var(58)*var(83)
    r56 = rct(57)*var(68)*var(83)
    r57 = rct(58)*var(64)*var(83)
    r58 = rct(59)*var(61)*var(83)
    r59 = rct(60)*var(61)*var(83)
    r60 = rct(61)*var(57)*var(83)
    r61 = rct(62)*var(57)*var(83)
    r62 = rct(63)*var(53)*var(83)
    r63 = rct(64)*var(66)*var(83)
    r64 = rct(65)*var(63)*var(83)
    r65 = rct(66)*var(62)*var(83)
    r66 = rct(67)*var(16)*var(85)
    r67 = rct(68)*var(75)*var(89)
    r68 = rct(69)*var(71)*var(89)
    r69 = rct(70)*var(87)*var(91)
    r70 = rct(71)*var(12)
    r71 = rct(72)*var(76)*var(87)
    r72 = rct(73)*var(44)
    r73 = rct(74)*var(83)*var(91)
    r74 = rct(75)*var(69)*var(83)
    r75 = rct(76)*var(76)*var(83)
    r76 = rct(77)*var(71)*var(85)
    r77 = rct(78)*var(51)*var(89)
    r78 = rct(79)*var(51)*var(89)
    r79 = rct(80)*var(55)*var(92)
    r80 = rct(81)*var(65)*var(92)
    r81 = rct(82)*var(74)*var(84)
    r82 = rct(83)*var(56)*var(84)
    r83 = rct(84)*var(67)*var(84)
    r84 = rct(85)*var(70)*var(84)
    r85 = rct(86)*var(78)*var(84)
    r86 = rct(87)*var(58)*var(84)
    r87 = rct(88)*var(68)*var(84)
    r88 = rct(89)*var(64)*var(84)
    r89 = rct(90)*var(61)*var(84)
    r90 = rct(91)*var(57)*var(84)
    r91 = rct(92)*var(53)*var(84)
    r92 = rct(93)*var(66)*var(84)
    r93 = rct(94)*var(63)*var(84)
    r94 = rct(95)*var(62)*var(84)
    r95 = rct(96)*var(81)*var(89)
    r96 = rct(97)*var(82)*var(92)
    r97 = rct(98)*var(81)*var(85)
    r98 = rct(99)*var(74)*var(92)
    r99 = rct(100)*var(56)*var(92)
    r100 = rct(101)*var(67)*var(92)
    r101 = rct(102)*var(70)*var(92)
    r102 = rct(103)*var(78)*var(92)
    r103 = rct(104)*var(58)*var(92)
    r104 = rct(105)*var(68)*var(92)
    r105 = rct(106)*var(64)*var(92)
    r106 = rct(107)*var(61)*var(92)
    r107 = rct(108)*var(57)*var(92)
    r108 = rct(109)*var(53)*var(92)
    r109 = rct(110)*var(66)*var(92)
    r110 = rct(111)*var(63)*var(92)
    r111 = rct(112)*var(62)*var(92)
    r112 = rct(113)*var(50)*var(89)
    r113 = rct(114)*var(82)**2
    r114 = rct(115)*var(82)**2
    r115 = rct(116)*var(82)*var(84)
    r116 = rct(117)*var(55)*var(84)
    r117 = rct(118)*var(65)*var(84)
    r118 = rct(119)*var(84)*var(88)
    r119 = rct(120)*var(84)*var(91)
    r120 = rct(121)*var(69)*var(84)
    r121 = rct(122)*var(76)*var(84)
    r122 = rct(123)*var(47)*var(89)
    r123 = rct(124)*var(47)*var(90)
    r124 = rct(125)*var(44)*var(89)
    r125 = rct(126)*var(44)*var(90)
    r126 = rct(127)*var(49)*var(89)
    r127 = rct(128)*var(47)*var(85)
    r128 = rct(129)*var(23)*var(89)
    r129 = rct(130)*var(73)*var(89)
    r130 = rct(131)*var(23)*var(85)
    r131 = rct(132)*var(73)*var(85)
    r132 = rct(133)*var(46)*var(89)
    r133 = rct(134)*var(77)*var(89)
    r134 = rct(135)*var(79)*var(89)
    r135 = rct(136)*var(60)*var(89)
    r136 = rct(137)*var(55)*var(88)
    r137 = rct(138)*var(65)*var(88)
    r138 = rct(139)*var(55)*var(88)
    r139 = rct(140)*var(65)*var(88)
    r140 = rct(141)*var(46)*var(90)
    r141 = rct(142)*var(77)*var(90)
    r142 = rct(143)*var(79)*var(90)
    r143 = rct(144)*var(46)*var(85)
    r144 = rct(145)*var(79)*var(85)
    r145 = rct(146)*var(79)*var(85)
    r146 = rct(147)*var(91)*var(92)
    r147 = rct(148)*var(69)*var(92)
    r148 = rct(149)*var(76)*var(92)
    r149 = rct(150)*var(91)*var(92)
    r150 = rct(151)*var(69)*var(92)
    r151 = rct(152)*var(76)*var(92)
    r152 = rct(153)*var(34)*var(89)
    r153 = rct(154)*var(33)*var(89)
    r154 = rct(155)*var(24)*var(89)
    r155 = rct(156)*var(26)*var(89)
    r156 = rct(157)*var(27)*var(89)
    r157 = rct(158)*var(25)*var(89)
    r158 = rct(159)*var(28)*var(89)
    r159 = rct(160)*var(32)*var(89)
    r160 = rct(161)*var(31)*var(89)
    r161 = rct(162)*var(45)*var(89)
    r162 = rct(163)*var(42)*var(89)
    r163 = rct(164)*var(41)*var(89)
    r164 = rct(165)*var(43)*var(89)
    r165 = rct(166)*var(38)*var(89)
    r166 = rct(167)*var(20)*var(89)
    r167 = rct(168)*var(21)*var(89)
    r168 = rct(169)*var(40)*var(85)
    r169 = rct(170)*var(59)*var(89)
    r170 = rct(171)*var(59)*var(90)
    r171 = rct(172)*var(88)**2
    r172 = rct(173)*var(88)*var(92)
    r173 = rct(174)*var(88)*var(92)
    r174 = rct(175)*var(74)*var(88)
    r175 = rct(176)*var(67)*var(88)
    r176 = rct(177)*var(70)*var(88)
    r177 = rct(178)*var(78)*var(88)
    r178 = rct(179)*var(58)*var(88)
    r179 = rct(180)*var(68)*var(88)
    r180 = rct(181)*var(64)*var(88)
    r181 = rct(182)*var(61)*var(88)
    r182 = rct(183)*var(57)*var(88)
    r183 = rct(184)*var(66)*var(88)
    r184 = rct(185)*var(56)*var(88)
    r185 = rct(186)*var(53)*var(88)
    r186 = rct(187)*var(63)*var(88)
    r187 = rct(188)*var(62)*var(88)
    r188 = rct(189)*var(74)*var(88)
    r189 = rct(190)*var(67)*var(88)
    r190 = rct(191)*var(70)*var(88)
    r191 = rct(192)*var(78)*var(88)
    r192 = rct(193)*var(58)*var(88)
    r193 = rct(194)*var(68)*var(88)
    r194 = rct(195)*var(61)*var(88)
    r195 = rct(196)*var(57)*var(88)
    r196 = rct(197)*var(56)*var(88)
    r197 = rct(198)*var(64)*var(88)
    r198 = rct(199)*var(53)*var(88)
    r199 = rct(200)*var(63)*var(88)
    r200 = rct(201)*var(62)*var(88)
    r201 = rct(202)*var(66)*var(88)
    r202 = rct(203)*var(82)*var(88)
    r203 = rct(204)*var(82)*var(88)
    r204 = rct(205)*var(88)*var(91)
    r205 = rct(206)*var(69)*var(88)
    r206 = rct(207)*var(76)*var(88)
    r207 = rct(208)*var(85)**2
    r208 = fix(3)*rct(209)
    r209 = fix(3)*rct(210)
    r210 = fix(3)*rct(211)
    r211 = fix(3)*rct(212)
    r212 = fix(3)*rct(213)
    r213 = fix(3)*rct(214)
    r214 = fix(3)*rct(215)
    r215 = fix(3)*rct(216)
    r216 = fix(3)*rct(217)
    r217 = fix(3)*rct(218)
    r218 = fix(3)*rct(219)
    r219 = fix(3)*rct(220)
    r220 = fix(3)*rct(221)
    r221 = fix(3)*rct(222)
    r222 = fix(3)*rct(223)
    r223 = fix(3)*rct(224)
    r224 = fix(3)*rct(225)
    r225 = rct(226)*var(87)
    r226 = rct(227)*var(90)
    r227 = rct(228)*var(15)
    r228 = rct(229)*var(54)
    r229 = rct(230)*var(72)
    r230 = rct(231)*var(18)
    r231 = rct(232)*var(11)
    r232 = rct(233)*var(44)
    r233 = rct(234)*var(12)
    r234 = rct(235)*var(75)
    r235 = rct(236)*var(77)
    r236 = rct(237)*var(79)
    r237 = rct(238)*var(84)
    r238 = rct(239)*var(87)
    r239 = rct(240)*var(85)
    r240 = rct(241)*var(18)
    r241 = rct(242)*var(29)*var(89)
    r242 = rct(243)*var(29)*var(89)
    r243 = rct(244)*var(29)*var(85)
    r244 = rct(245)*var(14)*var(89)
    r245 = rct(246)*var(86)*var(90)
    r246 = rct(247)*var(52)*var(84)
    r247 = rct(248)*var(84)*var(86)
    r248 = rct(249)*var(37)*var(89)
    r249 = rct(250)*var(52)**2
    r250 = rct(251)*var(52)**2
    r251 = rct(252)*var(52)*var(83)
    r252 = rct(253)*var(39)*var(86)
    r253 = rct(254)*var(36)*var(89)
    r254 = rct(255)*var(52)*var(89)
    r255 = rct(256)*var(85)*var(86)
    r256 = rct(257)*var(72)*var(86)
    r257 = rct(258)*var(80)*var(86)
    r258 = rct(259)*var(51)*var(86)
    r259 = rct(260)*var(40)*var(86)
    r260 = rct(261)*var(30)*var(86)
    r261 = rct(262)*var(86)*var(87)
    r262 = rct(263)*var(52)*var(87)
    r263 = rct(264)*var(8)*var(89)
    r264 = rct(265)*var(89)*var(9)
    r265 = rct(266)*var(10)*var(89)
    r266 = rct(267)*var(39)
    r267 = rct(268)*var(35)
    r268 = rct(269)*var(37)
    r269 = rct(270)*var(35)
    r270 = rct(271)*var(37)
    r271 = rct(272)*var(35)
    r272 = rct(273)*var(37)
    r273 = rct(274)*var(39)
    r274 = rct(275)*var(36)
    r275 = rct(276)*var(90)
    r276 = rct(277)*var(87)
    r277 = rct(278)*var(11)
    r278 = rct(279)*var(22)
    r279 = rct(280)*var(72)
    r280 = rct(281)*var(72)
    r281 = rct(282)*var(54)
    r282 = rct(283)*var(17)
    r283 = rct(284)*var(19)
    r284 = rct(285)*var(85)
    r285 = rct(286)*var(85)
    r286 = rct(287)*var(18)
    r287 = rct(288)*var(18)
    r288 = rct(289)*var(19)
    r289 = rct(290)*var(80)
    r290 = rct(291)*var(80)
    r291 = rct(292)*var(15)
    r292 = rct(293)*var(71)
    r293 = rct(294)*var(51)
    r294 = rct(295)*var(51)
    r295 = rct(296)*var(81)
    r296 = rct(297)*var(49)
    r297 = rct(298)*var(23)
    r298 = rct(299)*var(23)
    r299 = rct(300)*var(73)
    r300 = rct(301)*var(73)
    r301 = rct(302)*var(77)
    r302 = rct(303)*var(77)
    r303 = rct(304)*var(77)
    r304 = rct(305)*var(79)
    r305 = rct(306)*var(79)
    r306 = rct(307)*var(60)
    r307 = rct(308)*var(34)
    r308 = rct(309)*var(33)
    r309 = rct(310)*var(24)
    r310 = rct(311)*var(26)
    r311 = rct(312)*var(27)
    r312 = rct(313)*var(25)
    r313 = rct(314)*var(32)
    r314 = rct(315)*var(28)
    r315 = rct(316)*var(31)
    r316 = rct(317)*var(45)
    r317 = rct(318)*var(42)
    r318 = rct(319)*var(41)
    r319 = rct(320)*var(43)
    r320 = rct(321)*var(38)
    r321 = rct(322)*var(20)
    r322 = rct(323)*var(75)
    r323 = rct(324)*var(21)
    r324 = rct(325)*var(36)
    r325 = rct(326)*var(52)
    r326 = rct(327)*var(35)
    r327 = rct(328)*var(39)
    r328 = rct(329)*var(39)
    r329 = rct(330)*var(13)
    r330 = rct(331)*var(8)

    vardot(1) = r244
    vardot(2) = 0.25*r242
    vardot(3) = r12 + 0.15*r140 + 0.16*r142 + r40
    vardot(4) = r225 + r226 + r227 + r228 + r229 + r230 + r231 + r232 + r233 &
        + r234 + r235 + r236 + r271 + r272 + r273 + r274
    vardot(5) = 0
    vardot(6) = 0
    vardot(7) = r132
    vardot(8) = r222 - r263 - r330
    vardot(9) = r223 - r264
    vardot(10) = -r265
    vardot(11) = r10 - r231 + 0.5*r237 - r277 + r6 - r8
    vardot(12) = -r233 + r69 - r70
    vardot(13) = r261 - r329
    vardot(14) = r241 + 0.75*r242 + r243 - r244
    vardot(15) = -r227 - r291 + r38 - r39
    vardot(16) = r211 - r47 - r66
    vardot(17) = r23 + 0.5*r238 - r24 - r282
    vardot(18) = -r230 - r240 - r286 - r287 + r31 - r32
    vardot(19) = r25 - r26 - r27 - r283 - r288
    vardot(20) = 0.7*r121 - r166 - r321
    vardot(21) = 0.41*r118 - r167 - r323
    vardot(22) = r15 - r18 - r19 - r278
    vardot(23) = -r128 - r130 - r297 - r298
    vardot(24) = r115 - r154 - r309
    vardot(25) = -r157 - r312 + r81
    vardot(26) = r116 - r155 - r310
    vardot(27) = -r156 - r311 + r92
    vardot(28) = 0.7*r119 - r158 - r314
    vardot(29) = -r241 - r242 - r243
    vardot(30) = r215 - r260 - r43 - r44
    vardot(31) = 0.71*r120 - r160 - r315
    vardot(32) = r117 - r159 - r313
    vardot(33) = -r153 - r308 + r94
    vardot(34) = -r152 - r307 + r93
    vardot(35) = r246 + r253 + r266 - r267 - r269 - r271 - r326
    vardot(36) = r224 + r250 + r252 - r253 + 0.5*r267 + 0.5*r268 + 0.5*r269 + &
        0.5*r270 - r274 - r324
    vardot(37) = r247 - r248 + r256 + r257 + r258 + r259 + r260 - r268 - r270 &
        - r272
    vardot(38) = -r165 - r320 + r90
    vardot(39) = -r252 + r262 - r266 - r273 - r327 - r328
    vardot(40) = -r168 + r216 - r259 - r41
    vardot(41) = -r163 - r318 + r88 + r91
    vardot(42) = -r162 - r317 + r87
    vardot(43) = -r164 - r319 + r89
    vardot(44) = -r124 - r125 - r232 + r71 - r72
    vardot(45) = -r161 - r316 + r85 + r86
    vardot(46) = -r132 - r140 - r143 + r212
    vardot(47) = -r122 - r123 - r127 + 0.07*r140 + r214 + r301
    vardot(48) = 0.33*r104 + 0.15*r107 - r12 + 0.42*r123 + 0.4*r126 + &
        2.0*r128 + r129 + 2.0*r130 + r131 + 0.05*r140 + 0.05*r141 + 0.2*r142 &
        + 0.4*r170 + 0.65*r179 + 0.83*r182 + r20 + r210 + r256 + r257 + r258 &
        + r259 + r260 + r279 + r280 + r289 + r290 + r292 + r294 + r296 + &
        1.5*r297 + 2.0*r298 + r299 + r300 + r301 + r302 + r305 + 0.67*r317 + &
        0.5*r320 + r35 + 0.05*r36 + 0.61*r56
    vardot(49) = 0.13*r104 + 0.5*r105 + 0.36*r106 - r126 + 0.28*r170 + &
        0.26*r179 + r180 + 0.72*r181 - r296 + 0.26*r317 + 0.7*r319 + 0.24*r56 &
        + 0.95*r57 + 0.72*r58
    vardot(50) = 0.25*r101 + 0.25*r102 + 0.25*r103 + 0.25*r104 + 0.25*r105 + &
        0.25*r106 + 0.25*r108 + 0.25*r109 + 0.25*r110 + 0.25*r111 - r112 + &
        0.25*r79 + 0.25*r80 + 0.25*r98 + 0.25*r99
    vardot(51) = 0.75*r109 + 0.5*r156 + 0.32*r174 + r183 + r201 + r213 - r258 &
        - r293 - r294 + r311 + 0.32*r322 + 0.32*r48 + r63 - r77 - r78 + &
        0.16*r98
    vardot(52) = r245 - r246 - 2.0*r249 - 2.0*r250 - r251 - r254 + r255 - &
        r262 - r325 + r328
    vardot(53) = -r108 + r144 - r185 - r198 - r62 - r91
    vardot(54) = 0.425*r110 + r130 + r131 + r145 + r168 + 0.85*r186 + r21 - &
        r22 + r221 - r228 + 0.5*r238 + r239 + 2.0*r240 + r243 + r266 - r281 + &
        r35 + r37 + r55 + 0.08*r56 + 0.05*r57 + r59 + r61 + 0.85*r64 + r66 + &
        r76 + r97
    vardot(55) = -r116 - r136 - r138 + 0.5*r155 + 0.05*r174 + r260 + &
        0.05*r322 + r44 - r45 + 0.05*r48 - r79 + 0.03*r98
    vardot(56) = -r184 - r196 - r50 + r67 - r82 - r99
    vardot(57) = -r107 + 0.43*r134 + r165 - r182 - r195 - r60 - r61 - r90
    vardot(58) = 0.07*r102 - r103 + 0.136*r177 - r178 - r192 - r54 - r55 - &
        r86
    vardot(59) = 0.06*r102 + 0.5*r103 + 0.509*r161 - r169 - r170 + 0.127*r177 &
        + r178 + 0.373*r316 + 0.34*r53 + r54
    vardot(60) = 0.2*r100 + 0.18*r104 + 0.5*r105 + r107 + 0.59*r124 - r135 + &
        0.65*r139 + 0.2*r170 + 0.36*r179 + r180 + 0.83*r182 - r306 + &
        0.36*r317 + r320 + 0.33*r56 + 0.95*r57 + r60 + 0.16*r80
    vardot(61) = -r106 + r133 + 0.5*r164 - r181 - r194 - r58 - r59 - r89
    vardot(62) = -r111 + r127 + r153 - r187 - r200 - r65 - r94
    vardot(63) = -r110 + r143 + r152 - r186 - r199 - r64 - r93
    vardot(64) = -r105 + 0.5*r163 - r180 - r197 - r57 - r88
    vardot(65) = -r117 + r122 - r137 - r139 + r159 - r46 - r80
    vardot(66) = -r109 + 0.5*r156 + 0.18*r174 - r183 - r201 + 0.18*r322 + r43 &
        + 0.18*r48 - r63 - r92 + 0.09*r98
    vardot(67) = -r100 - r175 - r189 + r258 - r51 + r77 + r78 - r83
    vardot(68) = -r104 + r162 + 0.44*r169 - r179 - r193 - r56 - r87
    vardot(69) = -r120 + 0.8*r126 - r147 - r150 + r160 - r205 - r74
    vardot(70) = -r101 - r176 - r190 - r52 - r84 + r95 + r97
    vardot(71) = 0.25*r105 + 0.25*r108 + 0.25*r110 + 0.25*r111 + r112 + r136 &
        + r138 + 0.35*r139 + 0.5*r155 + 0.5*r157 + 0.5*r163 + 0.5*r164 + &
        0.13*r174 + 0.57*r184 + r196 + r197 + r198 + r199 + r200 - r292 + &
        r307 + r308 + r310 + r312 + r318 + 0.13*r322 + r45 + 0.13*r48 + &
        0.57*r50 - r68 - r76 + 0.75*r79 + 0.09*r80 + 0.07*r98 + 0.54*r99
    vardot(72) = 0.5*r100 + 0.75*r101 + 1.1*r102 + 1.13*r103 + 0.95*r104 + &
        0.75*r105 + 0.89*r106 + 0.85*r107 + 1.25*r108 + 0.75*r109 + 0.83*r110 &
        + 1.25*r111 + 0.29*r120 + 0.535*r123 + 2.23*r124 + 0.6*r125 + r137 + &
        r14 + 0.9*r140 + 0.8*r141 + 0.7*r142 + r146 + 2.0*r147 + 2.0*r148 + &
        r149 + r150 + r151 + r16 + 0.5*r167 + 2.0*r17 + 0.12*r170 + r172 + &
        r173 + 0.2*r175 + 0.69*r177 + 0.75*r178 + 0.4*r179 + 0.28*r181 + &
        0.17*r182 + 0.39*r184 + r185 + 0.15*r186 + r187 + r19 - r20 + r205 + &
        r206 + r219 - r229 + r241 + r243 - r256 + r278 - r279 - r280 + r296 + &
        0.5*r297 + r302 + r305 + r306 + r313 + r315 + 0.627*r316 + 0.3*r319 + &
        0.5*r320 + r321 + r33 - r35 + 0.05*r36 + r4 + r46 + 0.39*r50 + &
        0.96*r51 + 0.56*r53 + 0.75*r54 + 0.35*r56 + 0.28*r58 + r60 + r62 + &
        0.15*r64 + r65 + r74 + r75 + 0.75*r79 + 1.25*r80 + 0.75*r96 + &
        0.75*r98 + 0.95*r99
    vardot(73) = 0.5*r100 + 0.29*r104 + 0.14*r106 + 0.5*r108 - r129 - r131 + &
        r135 + 0.82*r141 + 0.8*r142 + 0.6*r170 + 0.8*r175 + 0.58*r179 + &
        0.28*r181 + 0.17*r182 + r185 - r299 - r300 + 0.58*r317 + 0.3*r319 + &
        0.53*r56 + 0.28*r58 + r62 + r84
    vardot(74) = 0.5*r157 - r174 + 0.3*r184 - r188 + r47 - r48 - r49 + &
        0.3*r50 + r66 - r81 - r98 + 0.15*r99
    vardot(75) = -r234 - r322 + r49 + 0.04*r51 + 0.07*r52 - r67 + r82
    vardot(76) = -r121 + 0.57*r134 + r145 - r148 - r151 + r166 + 0.41*r169 - &
        r206 + r303 + r304 - r71 + r72 - r75
    vardot(77) = 0.2*r102 + 0.03*r110 - r133 + 0.159*r140 - r141 + 0.402*r177 &
        + 0.05*r186 - r235 - r301 - r302 - r303 + 0.368*r316 + 0.34*r53 + &
        0.05*r64
    vardot(78) = -r102 + r132 + 0.491*r161 - r177 - r191 - r53 - r85
    vardot(79) = 0.14*r102 + 0.05*r110 - r134 + 0.387*r140 - r142 - r144 - &
        r145 + 0.288*r177 + 0.1*r186 - r236 - r304 - r305 + 0.259*r316 + &
        0.22*r53 + 0.1*r64
    vardot(80) = 0.5*r101 + 0.5*r111 + 2.0*r113 + r114 + 0.5*r123 + r137 + &
        0.04*r141 + 0.5*r154 + 0.5*r158 + 0.32*r174 + r176 + 0.75*r184 + r187 &
        + r202 + r203 + r218 - r257 - r289 - r290 + r300 + r309 + r313 + r314 &
        + 0.32*r322 - r36 - r37 + r42 + r46 + 0.32*r48 + 0.75*r50 + 0.93*r52 &
        + r65 + 0.5*r80 + 0.75*r96 + 0.16*r98 + 0.38*r99
    vardot(81) = 0.25*r101 + 0.25*r102 + 0.25*r103 + 0.25*r104 + 0.25*r106 + &
        0.19*r174 + r188 + r189 + r190 + r191 + r192 + r193 + r194 + r195 + &
        r217 - r295 + 0.19*r322 + 0.19*r48 - r95 - r97 + 0.35*r98
    vardot(82) = -2.0*r113 - 2.0*r114 - r115 + r146 + 0.5*r154 + r168 + &
        0.32*r174 - r202 - r203 + r204 + r259 + r292 + 0.85*r295 + 0.32*r322 &
        + r41 - r42 + 0.32*r48 + r73 - r96 + 0.16*r98
    vardot(83) = -r0 - r14 + r208 - r23 - r251 + r276 + r282 + r285 + r287 - &
        r29 + r34 - r40 - r42 - r45 - r46 - r48 - r49 - r50 - r51 - r52 - r53 &
        - r54 - r55 - r56 - r57 - r58 - r59 - r60 - r61 - r62 - r63 - r64 - &
        r65 - r73 - r74 - r75 - r9
    vardot(84) = r1 - 2.0*r10 + 0.3*r100 + 0.5*r101 + 0.92*r102 + r103 + r104 &
        + 0.5*r105 + 0.64*r106 + 1.15*r107 + 0.5*r108 + r109 + r11 + &
        0.45*r110 + 0.5*r111 + r112 + 2.0*r113 - r115 - r116 - r117 - r118 - &
        r119 + r12 - r120 - r121 + 0.3*r123 + 2.0*r124 + r125 + 0.2*r126 + &
        r128 + r130 + r135 + r136 + r137 + r14 + 0.06*r140 + 0.06*r141 + &
        0.275*r142 + r146 + 2.0*r147 + r148 - r15 + 0.15*r169 + 2.0*r17 + &
        r172 + 0.27*r174 + 0.8*r175 + 0.864*r177 + r178 + r179 + 0.28*r181 + &
        r182 + r183 + 0.8*r186 - r2 + r20 + r202 + r205 - r237 + r244 - r246 &
        - r247 - r25 + r254 + r256 + r26 + r278 + 2.0*r279 - r28 + r288 + &
        r289 + r292 + 2.0*r296 + 2.0*r298 + r299 + r30 + r302 + r304 + r305 + &
        r306 + r307 + r308 + r309 + r310 + r311 + r312 + r313 + r314 + r315 + &
        r316 + r317 + r318 + 0.3*r319 + r320 + 0.27*r322 + r33 + r35 + &
        0.05*r36 + r4 + r42 + r45 + r46 + 0.27*r48 + 0.9*r53 + r54 + 0.92*r56 &
        + 0.05*r57 + 0.28*r58 + r60 + r63 + 0.8*r64 - r7 + r74 + r79 + r8 + &
        r80 - r81 - r82 - r83 - r84 - r85 - r86 - r87 - r88 - r89 - r9 - r90 &
        - r91 - r92 - r93 - r94 + r96 + 0.64*r98 + 0.5*r99
    vardot(85) = -r127 - r130 - r131 - r143 - r144 - r145 - r168 - 2.0*r207 + &
        r22 - r239 - r243 + r252 - r255 - r28 + r283 - r284 - r285 + r286 + &
        r287 - r29 + 0.4*r291 + r3 - r30 - r31 + r32 + r327 - r34 - r35 - r37 &
        - r66 - r76 - r97
    vardot(86) = -r245 - r247 + r248 + 2.0*r249 + r251 - r252 + r253 + r254 - &
        r255 - r256 - r257 - r258 - r259 - r260 - r261 + 3.0*r263 + 2.0*r264 &
        + r265 + 2.0*r324 + r325 + r326 + r327 + r329 + 3.0*r330
    vardot(87) = r0 + r105 + r108 + 0.575*r110 + r111 + r124 + r125 + r14 + &
        0.5*r163 + r180 + r184 + r185 + 0.15*r186 + r187 + r196 + r197 + r198 &
        + r199 + r200 + 2.0*r207 + r209 - r21 - r225 - r238 + r24 - r25 + &
        r251 + r255 + r26 - r261 - r262 + r27 - r276 + r28 + r281 + r284 + &
        r286 + r288 + 2.0*r29 + 0.6*r291 - r3 + r30 + r307 + r308 - r31 + &
        r318 + r32 + r322 + r328 + r329 - r38 + r39 + r40 + r42 + r45 + r46 + &
        r48 + 2.0*r50 + 0.96*r51 + 0.93*r52 + 0.9*r53 + r54 + 0.92*r56 + &
        1.95*r57 + r58 + r60 + 2.0*r62 + r63 + 1.15*r64 + 2.0*r65 - r69 + r70 &
        - r71 + r72 + r73 + r74 + r75 + r9 + r99
    vardot(88) = 0.3*r100 + 0.5*r101 + 0.36*r106 - r118 + r129 + r131 - r136 &
        - r137 - r138 - r139 + r148 + 0.5*r167 - 2.0*r171 - r172 - r173 - &
        r174 - 0.8*r175 - r177 - r178 - r179 - r180 - 0.28*r181 - r182 - r183 &
        - r184 - r185 - r186 - r187 - r188 - r189 - r190 - r191 - r192 - r193 &
        - r194 - r195 - r196 - r197 - r198 - r199 - r200 - r201 - r202 - r203 &
        - r204 - r205 + r257 + 0.6*r291 + r293 + 0.85*r295 + r299 + r302 + &
        r305 + r306 + 0.7*r319 + r321 + 0.95*r36 + r37 - r38 + r39 - r40 + &
        0.96*r51 + 0.93*r52 + 0.72*r58 + r75 + r83
    vardot(89) = -r1 - r11 - r112 + 0.44*r118 - r12 - r122 + 0.135*r123 - &
        r124 - r126 - r128 - r129 - r13 - r132 - r133 - r134 - r135 + &
        0.27*r140 + 0.08*r141 + 0.215*r142 - r152 - r153 - 0.5*r154 - &
        0.5*r155 - 0.5*r156 - 0.5*r157 - 0.5*r158 - r159 - r160 - 0.491*r161 &
        - r162 - 0.5*r163 - 0.5*r164 - r165 - r166 - 0.5*r167 - r169 + &
        0.1*r170 - r18 + r2 - r20 - r21 - r22 - r23 - r24 - r241 - r242 - &
        r244 - r248 - r253 - r254 - r263 - r264 - r265 - r27 + 2.0*r275 + &
        2.0*r277 + r278 + r28 + r281 + r282 + r283 - r30 + r307 + r308 + r309 &
        + r310 + r311 + r312 + r313 + r314 + r315 + r316 + r317 + r318 + r319 &
        + r320 + r321 + r323 + r326 - r33 - r36 - r41 - r43 - r44 - r47 - &
        2.0*r5 - 2.0*r6 - r67 - r68 - r7 - r77 - r78 - r8 + r9 - r95
    vardot(90) = -r0 - r1 + 0.15*r118 + 0.3*r119 + 0.29*r120 + 0.3*r121 - &
        r123 - r125 - 0.9*r140 - 0.8*r141 - 0.8*r142 - 0.7*r170 - r2 + r220 - &
        r226 - r245 - r275 + r276 + r284 + r287 - r3 + r325 - r4 + r5
    vardot(91) = -r119 - r146 - r149 + 0.5*r158 - r204 + 0.15*r295 + r68 - &
        r69 + r70 - r73 + r76
    vardot(92) = -r100 - r101 - r102 - r103 - r104 - r105 - r106 - r107 - &
        r108 - r109 - r110 - r111 + 0.44*r118 + 0.305*r123 + r13 + r136 + &
        r137 - r14 - r146 - r147 - r148 - r149 - r15 - r150 - r151 - 2.0*r16 &
        - 2.0*r17 + 2.0*r171 - r173 + 1.18*r174 + r175 + r176 + r177 + r178 + &
        r179 + r18 + r180 + r181 + r182 + r183 + r184 + r185 + r186 + r187 + &
        r202 + r204 + r205 + r206 + r241 + r242 + r243 + r289 + 0.4*r291 + &
        r293 + 2.0*r294 + 0.15*r295 + r303 + 0.18*r322 + r323 - r4 + r40 + &
        0.18*r48 - r79 - r80 + r83 + r84 - r96 - 0.91*r98 - r99

  END SUBROUTINE Fun



END MODULE gckpp_function
!------------------------- END gckpp_function.f90 END ------------------------
