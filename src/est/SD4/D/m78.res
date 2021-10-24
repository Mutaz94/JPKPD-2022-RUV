Sun Oct 24 04:24:44 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/SD4/D/dat78.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2

$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr
$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m78.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1634.60401528236        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2480E+02 -9.0199E+01 -4.8500E+01 -7.7433E+01  2.9304E+01  3.1646E+01 -3.6463E+01  1.3013E+01 -4.8807E+01  7.9255E+00
             9.7504E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1658.07153531831        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.9522E-01  1.1356E+00  1.2457E+00  1.0105E+00  1.1309E+00  1.0165E+00  1.2247E+00  9.2760E-01  1.2605E+00  9.7239E-01
             9.6965E-01
 PARAMETER:  9.5212E-02  2.2712E-01  3.1968E-01  1.1044E-01  2.2300E-01  1.1635E-01  3.0268E-01  2.4840E-02  3.3151E-01  7.1997E-02
             6.9184E-02
 GRADIENT:   1.4408E+01 -1.6870E+01 -3.1060E+00 -1.0030E+01 -9.3767E+00  2.9342E+00 -4.1306E+00  8.4191E-01  1.3577E+01 -1.2918E+01
            -8.0072E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1659.52218251856        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.7402E-01  1.2770E+00  1.4584E+00  9.5826E-01  1.2975E+00  9.9186E-01  1.2950E+00  9.1703E-01  1.1883E+00  1.2058E+00
             1.0050E+00
 PARAMETER:  7.3677E-02  3.4455E-01  4.7732E-01  5.7361E-02  3.6045E-01  9.1831E-02  3.5852E-01  1.3384E-02  2.7250E-01  2.8712E-01
             1.0497E-01
 GRADIENT:  -3.4228E+01  1.4340E+01  9.6068E+00  8.7661E+00  2.0940E+00 -7.3445E+00  4.0694E+00 -4.1073E+00 -9.2894E-01  1.6646E+00
             3.4927E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1660.89810842748        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.9086E-01  1.3118E+00  1.3961E+00  9.2594E-01  1.2828E+00  1.0081E+00  1.2046E+00  1.1601E+00  1.2357E+00  1.1609E+00
             9.8494E-01
 PARAMETER:  9.0816E-02  3.7140E-01  4.3367E-01  2.3053E-02  3.4907E-01  1.0807E-01  2.8618E-01  2.4849E-01  3.1166E-01  2.4918E-01
             8.4822E-02
 GRADIENT:   3.6230E+00  3.1435E+00  3.2367E-01  4.3977E+00 -1.4025E+00 -2.3786E-01  1.3857E-01  2.0730E-01 -2.7879E-02 -9.6672E-02
            -3.8502E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1660.96713688737        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      704
 NPARAMETR:  9.9190E-01  1.4894E+00  1.2218E+00  8.0862E-01  1.3082E+00  1.0096E+00  1.1030E+00  1.0443E+00  1.3429E+00  1.1676E+00
             9.8714E-01
 PARAMETER:  9.1868E-02  4.9835E-01  3.0032E-01 -1.1243E-01  3.6867E-01  1.0956E-01  1.9802E-01  1.4330E-01  3.9485E-01  2.5495E-01
             8.7060E-02
 GRADIENT:   3.7690E+00  6.1293E+00  2.6213E+00  4.5716E+00 -4.5307E+00  6.1448E-02 -1.8439E+00 -7.4525E-01 -6.5527E-01  6.7827E-02
            -3.4782E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1660.96858342748        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      881
 NPARAMETR:  9.9153E-01  1.5371E+00  1.1648E+00  7.7855E-01  1.3125E+00  1.0096E+00  1.0909E+00  9.9648E-01  1.3658E+00  1.1670E+00
             9.8709E-01
 PARAMETER:  9.1496E-02  5.2987E-01  2.5251E-01 -1.5032E-01  3.7192E-01  1.0958E-01  1.8702E-01  9.6473E-02  4.1175E-01  2.5445E-01
             8.7006E-02
 GRADIENT:   2.4050E+00  8.5410E+00  3.0495E+00  5.6602E+00 -5.2195E+00 -1.8838E-02 -1.9049E+00 -7.9967E-01 -7.2858E-01  5.3869E-02
            -4.7780E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1660.96935550190        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  9.9120E-01  1.5683E+00  1.1235E+00  7.5865E-01  1.3143E+00  1.0096E+00  1.0848E+00  9.5677E-01  1.3801E+00  1.1658E+00
             9.8699E-01
 PARAMETER:  9.1156E-02  5.5000E-01  2.1648E-01 -1.7621E-01  3.7333E-01  1.0960E-01  1.8140E-01  5.5807E-02  4.2212E-01  2.5341E-01
             8.6901E-02
 GRADIENT:   1.2936E+00  9.9954E+00  3.2026E+00  6.2857E+00 -5.5122E+00 -7.6542E-02 -1.8473E+00 -7.8476E-01 -7.5539E-01  4.3592E-02
            -5.4445E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1660.96998863299        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1240
 NPARAMETR:  9.9094E-01  1.5902E+00  1.0934E+00  7.4457E-01  1.3154E+00  1.0097E+00  1.0808E+00  9.2522E-01  1.3901E+00  1.1647E+00
             9.8690E-01
 PARAMETER:  9.0902E-02  5.6387E-01  1.8932E-01 -1.9495E-01  3.7412E-01  1.0962E-01  1.7773E-01  2.2273E-02  4.2941E-01  2.5246E-01
             8.6813E-02
 GRADIENT:   4.7698E-01  1.0876E+01  3.2404E+00  6.6534E+00 -5.6223E+00 -1.1483E-01 -1.7708E+00 -7.5339E-01 -7.6446E-01  3.6503E-02
            -5.7933E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1660.97031660662        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1418
 NPARAMETR:  9.9074E-01  1.6073E+00  1.0697E+00  7.3349E-01  1.3162E+00  1.0097E+00  1.0777E+00  8.9901E-01  1.3983E+00  1.1638E+00
             9.8683E-01
 PARAMETER:  9.0701E-02  5.7457E-01  1.6740E-01 -2.0994E-01  3.7473E-01  1.0964E-01  1.7480E-01 -6.4659E-03  4.3529E-01  2.5171E-01
             8.6746E-02
 GRADIENT:  -1.6423E-01  1.1475E+01  3.2357E+00  6.8963E+00 -5.6579E+00 -1.4362E-01 -1.6968E+00 -7.2029E-01 -7.6662E-01  3.1043E-02
            -6.0041E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1660.97046397127        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1598
 NPARAMETR:  9.9057E-01  1.6226E+00  1.0487E+00  7.2355E-01  1.3170E+00  1.0097E+00  1.0747E+00  8.7475E-01  1.4060E+00  1.1631E+00
             9.8678E-01
 PARAMETER:  9.0527E-02  5.8401E-01  1.4752E-01 -2.2359E-01  3.7534E-01  1.0966E-01  1.7207E-01 -3.3819E-02  4.4076E-01  2.5107E-01
             8.6691E-02
 GRADIENT:  -7.2482E-01  1.1937E+01  3.2100E+00  7.0771E+00 -5.6568E+00 -1.6739E-01 -1.6232E+00 -6.8605E-01 -7.6491E-01  2.6370E-02
            -6.1479E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1660.97056120298        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1778
 NPARAMETR:  9.9043E-01  1.6359E+00  1.0304E+00  7.1480E-01  1.3177E+00  1.0097E+00  1.0720E+00  8.5287E-01  1.4130E+00  1.1625E+00
             9.8674E-01
 PARAMETER:  9.0384E-02  5.9220E-01  1.2992E-01 -2.3575E-01  3.7592E-01  1.0968E-01  1.6957E-01 -5.9143E-02  4.4574E-01  2.5055E-01
             8.6649E-02
 GRADIENT:  -1.1928E+00  1.2279E+01  3.1726E+00  7.2051E+00 -5.6321E+00 -1.8665E-01 -1.5551E+00 -6.5330E-01 -7.6066E-01  2.2491E-02
            -6.2398E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1660.97059431061        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1959
 NPARAMETR:  9.9032E-01  1.6469E+00  1.0154E+00  7.0754E-01  1.3185E+00  1.0097E+00  1.0697E+00  8.3441E-01  1.4191E+00  1.1620E+00
             9.8671E-01
 PARAMETER:  9.0273E-02  5.9892E-01  1.1528E-01 -2.4596E-01  3.7646E-01  1.0970E-01  1.6740E-01 -8.1025E-02  4.5001E-01  2.5016E-01
             8.6619E-02
 GRADIENT:  -1.5637E+00  1.2519E+01  3.1332E+00  7.2911E+00 -5.5975E+00 -2.0206E-01 -1.4972E+00 -6.2487E-01 -7.5539E-01  1.9496E-02
            -6.2938E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1660.97064364450        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2141
 NPARAMETR:  9.9021E-01  1.6587E+00  9.9961E-01  6.9979E-01  1.3193E+00  1.0098E+00  1.0672E+00  8.1434E-01  1.4258E+00  1.1616E+00
             9.8668E-01
 PARAMETER:  9.0166E-02  6.0601E-01  9.9614E-02 -2.5697E-01  3.7709E-01  1.0972E-01  1.6502E-01 -1.0537E-01  4.5470E-01  2.4979E-01
             8.6594E-02
 GRADIENT:  -1.9311E+00  1.2731E+01  3.0839E+00  7.3611E+00 -5.5494E+00 -2.1554E-01 -1.4350E+00 -5.9359E-01 -7.4807E-01  1.6669E-02
            -6.3290E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1661.11192282295        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2328
 NPARAMETR:  9.9141E-01  1.6508E+00  9.8605E-01  6.8918E-01  1.3232E+00  1.0103E+00  1.0747E+00  8.4678E-01  1.4284E+00  1.1605E+00
             9.8619E-01
 PARAMETER:  9.1376E-02  6.0128E-01  8.5947E-02 -2.7226E-01  3.8005E-01  1.1023E-01  1.7201E-01 -6.6318E-02  4.5652E-01  2.4884E-01
             8.6091E-02
 GRADIENT:   7.3550E-01 -5.3114E+00  2.7187E-01 -9.2285E-01 -7.4632E-01  3.0586E-02 -2.6141E-03 -2.5196E-03  9.4365E-02  1.3317E-01
             2.1052E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1661.11315614287        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2506
 NPARAMETR:  9.9183E-01  1.6524E+00  9.8440E-01  6.8922E-01  1.3238E+00  1.0105E+00  1.0751E+00  8.4608E-01  1.4283E+00  1.1604E+00
             9.8577E-01
 PARAMETER:  9.1797E-02  6.0222E-01  8.4275E-02 -2.7219E-01  3.8050E-01  1.1048E-01  1.7237E-01 -6.7136E-02  4.5650E-01  2.4875E-01
             8.5670E-02
 GRADIENT:   1.6237E+00 -4.4650E+00  4.0776E-02 -6.9850E-02 -1.6687E-01  1.2002E-01  5.5551E-02  1.3074E-02  1.3953E-01  5.3778E-02
             3.8526E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1661.11327864745        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2686             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9183E-01  1.6528E+00  9.8364E-01  6.8895E-01  1.3238E+00  1.0105E+00  1.0750E+00  8.4282E-01  1.4281E+00  1.1601E+00
             9.8570E-01
 PARAMETER:  9.1801E-02  6.0245E-01  8.3509E-02 -2.7259E-01  3.8048E-01  1.1048E-01  1.7231E-01 -7.0999E-02  4.5637E-01  2.4854E-01
             8.5600E-02
 GRADIENT:   4.2444E+02  5.0419E+02  6.5063E-01  8.2451E+01  1.8346E+01  4.1600E+01  1.3956E+01  2.7232E-02  2.1377E+01  2.2032E+00
             7.8386E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1661.11334988544        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2865
 NPARAMETR:  9.9184E-01  1.6530E+00  9.8311E-01  6.8879E-01  1.3237E+00  1.0105E+00  1.0750E+00  8.4198E-01  1.4282E+00  1.1601E+00
             9.8570E-01
 PARAMETER:  9.1803E-02  6.0258E-01  8.2963E-02 -2.7282E-01  3.8044E-01  1.1048E-01  1.7229E-01 -7.1996E-02  4.5641E-01  2.4847E-01
             8.5596E-02
 GRADIENT:   1.6281E+00 -4.4392E+00  9.4479E-02 -9.7415E-02 -1.1681E-01  1.1765E-01  3.0421E-02 -6.0065E-03  9.0462E-02  9.2793E-03
            -1.4932E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1661.11344636683        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     3045             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9184E-01  1.6533E+00  9.8196E-01  6.8855E-01  1.3237E+00  1.0105E+00  1.0749E+00  8.4214E-01  1.4285E+00  1.1600E+00
             9.8572E-01
 PARAMETER:  9.1805E-02  6.0279E-01  8.1794E-02 -2.7317E-01  3.8041E-01  1.1049E-01  1.7226E-01 -7.1809E-02  4.5660E-01  2.4840E-01
             8.5619E-02
 GRADIENT:   4.2431E+02  5.0466E+02  5.0497E-01  8.2615E+01  1.8507E+01  4.1575E+01  1.3970E+01  5.1730E-02  2.1408E+01  2.2209E+00
             8.2789E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1661.11351936003        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     3224
 NPARAMETR:  9.9184E-01  1.6536E+00  9.8159E-01  6.8838E-01  1.3236E+00  1.0105E+00  1.0749E+00  8.4113E-01  1.4286E+00  1.1599E+00
             9.8573E-01
 PARAMETER:  9.1807E-02  6.0293E-01  8.1417E-02 -2.7341E-01  3.8036E-01  1.1049E-01  1.7225E-01 -7.3005E-02  4.5668E-01  2.4834E-01
             8.5623E-02
 GRADIENT:   1.5901E+00 -4.5464E+00 -1.5064E-02 -2.7865E-02 -1.4351E-02  1.1003E-01  4.5099E-02  1.3118E-02  1.3419E-01  2.5918E-02
             2.4996E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1661.11361673313        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3406             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9184E-01  1.6539E+00  9.8112E-01  6.8813E-01  1.3235E+00  1.0105E+00  1.0748E+00  8.3887E-01  1.4285E+00  1.1598E+00
             9.8569E-01
 PARAMETER:  9.1810E-02  6.0314E-01  8.0942E-02 -2.7378E-01  3.8028E-01  1.1049E-01  1.7216E-01 -7.5700E-02  4.5665E-01  2.4823E-01
             8.5585E-02
 GRADIENT:   4.2430E+02  5.0539E+02  6.2264E-01  8.2645E+01  1.8357E+01  4.1570E+01  1.3964E+01  3.0529E-02  2.1352E+01  2.2002E+00
             7.8903E-01

0ITERATION NO.:   97    OBJECTIVE VALUE:  -1661.11361673313        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     3465
 NPARAMETR:  9.9184E-01  1.6539E+00  9.8112E-01  6.8813E-01  1.3235E+00  1.0105E+00  1.0748E+00  8.3887E-01  1.4285E+00  1.1598E+00
             9.8569E-01
 PARAMETER:  9.1810E-02  6.0314E-01  8.0942E-02 -2.7378E-01  3.8028E-01  1.1049E-01  1.7216E-01 -7.5700E-02  4.5665E-01  2.4823E-01
             8.5585E-02
 GRADIENT:  -5.6205E-04 -9.1571E-02  7.5943E-02  6.5722E-02 -1.7875E-02 -4.2076E-04  6.2570E-04 -1.7441E-03 -1.9937E-02  1.2966E-03
            -1.1881E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     3465
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.6892E-05 -2.2747E-02 -2.2381E-02  1.8770E-02 -4.2723E-02
 SE:             2.9799E-02  2.3153E-02  7.9589E-03  2.1958E-02  2.2021E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9848E-01  3.2586E-01  4.9221E-03  3.9265E-01  5.2367E-02

 ETASHRINKSD(%)  1.6933E-01  2.2435E+01  7.3337E+01  2.6438E+01  2.6227E+01
 ETASHRINKVR(%)  3.3837E-01  3.9837E+01  9.2891E+01  4.5886E+01  4.5576E+01
 EBVSHRINKSD(%)  4.3796E-01  2.0828E+01  7.7666E+01  2.9833E+01  2.3100E+01
 EBVSHRINKVR(%)  8.7401E-01  3.7318E+01  9.5012E+01  5.0765E+01  4.0864E+01
 RELATIVEINF(%)  9.8904E+01  3.5775E+00  6.8780E-01  2.7718E+00  1.7409E+01
 EPSSHRINKSD(%)  4.4037E+01
 EPSSHRINKVR(%)  6.8682E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1661.1136167331281     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -925.96279016938990     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.13
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1661.114       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.92E-01  1.65E+00  9.81E-01  6.88E-01  1.32E+00  1.01E+00  1.07E+00  8.39E-01  1.43E+00  1.16E+00  9.86E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
 Elapsed finaloutput time in seconds:     0.00
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,      111.848
Stop Time:
Sun Oct 24 04:25:03 CDT 2021
