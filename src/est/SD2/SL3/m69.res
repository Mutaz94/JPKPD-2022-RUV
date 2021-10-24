Sat Oct 23 19:27:24 CDT 2021
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
$DATA ../../../../data/SD2/SL3/dat69.csv ignore=@
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
Current Date:       23 OCT 2021
Days until program expires : 176
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
 NO. OF DATA RECS IN DATA SET:      795
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

 TOT. NO. OF OBS RECS:      695
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
 RAW OUTPUT FILE (FILE): m69.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1657.36337117510        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4513E+02 -1.4632E+01 -1.2957E+02  3.1148E+02  2.4981E+02  4.2872E+01 -5.1279E+01 -3.3073E+01 -8.7587E+01 -5.1556E+01
            -2.3948E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2377.57690813853        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0334E+00  1.2434E+00  1.3641E+00  8.6803E-01  1.1154E+00  9.3859E-01  1.0675E+00  9.2219E-01  1.1928E+00  1.2948E+00
             1.8827E+00
 PARAMETER:  1.3285E-01  3.1782E-01  4.1048E-01 -4.1528E-02  2.0917E-01  3.6625E-02  1.6528E-01  1.8991E-02  2.7628E-01  3.5837E-01
             7.3270E-01
 GRADIENT:   1.8117E+02  9.9949E+01 -1.8577E+01  9.5397E+01  4.0028E+01 -2.3874E+01  1.3978E+01 -3.7907E+00  9.7764E+00  6.7940E+00
            -1.3831E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2379.97122241100        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      203
 NPARAMETR:  1.0255E+00  1.2585E+00  1.7159E+00  8.6374E-01  1.2202E+00  9.7481E-01  8.8570E-01  9.5750E-01  1.1526E+00  1.4580E+00
             1.8851E+00
 PARAMETER:  1.2519E-01  3.2993E-01  6.3995E-01 -4.6485E-02  2.9905E-01  7.4490E-02 -2.1382E-02  5.6566E-02  2.4204E-01  4.7707E-01
             7.3398E-01
 GRADIENT:  -3.0893E+01  2.1904E+01 -1.4892E+01  7.2677E+01  3.2862E+01 -2.3222E+01 -6.8507E+00 -8.8189E+00 -8.3077E+00  9.5249E+00
            -1.4709E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2395.87645234757        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      379
 NPARAMETR:  1.0481E+00  1.3278E+00  2.3071E+00  7.8274E-01  1.2921E+00  1.0420E+00  8.1747E-01  2.1356E+00  1.2162E+00  1.4456E+00
             2.0167E+00
 PARAMETER:  1.4701E-01  3.8354E-01  9.3599E-01 -1.4496E-01  3.5624E-01  1.4112E-01 -1.0154E-01  8.5872E-01  2.9577E-01  4.6851E-01
             8.0144E-01
 GRADIENT:   1.5633E+01 -1.1484E+01 -6.2676E+00  2.6502E+00  1.1071E+01  4.8222E+00 -2.2424E+00 -2.2893E-01 -2.7592E+00 -3.7885E-01
            -1.4045E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2397.01682707893        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      554
 NPARAMETR:  1.0402E+00  1.4598E+00  2.9780E+00  7.0029E-01  1.3728E+00  1.0288E+00  7.8084E-01  2.6185E+00  1.3397E+00  1.5422E+00
             2.0300E+00
 PARAMETER:  1.3944E-01  4.7829E-01  1.1913E+00 -2.5626E-01  4.1682E-01  1.2842E-01 -1.4739E-01  1.0626E+00  3.9248E-01  5.3320E-01
             8.0801E-01
 GRADIENT:  -4.8175E-01 -4.3901E-01  6.6884E-01  6.6241E-01 -1.4692E+00 -4.6827E-02 -3.8051E-01 -1.7708E-01 -4.7304E-03  2.1519E-02
            -6.0157E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2397.07115464664        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      733
 NPARAMETR:  1.0427E+00  1.6137E+00  2.5959E+00  6.0109E-01  1.3965E+00  1.0298E+00  7.9749E-01  2.6656E+00  1.4209E+00  1.5639E+00
             2.0273E+00
 PARAMETER:  1.4186E-01  5.7851E-01  1.0539E+00 -4.0902E-01  4.3394E-01  1.2939E-01 -1.2629E-01  1.0804E+00  4.5127E-01  5.4716E-01
             8.0670E-01
 GRADIENT:   4.0106E+00  7.9763E+00  1.5257E+00  2.8670E+00 -3.4452E+00  1.6286E-01  2.3306E-01 -2.1928E-02 -1.0066E+00 -7.0347E-01
            -2.1392E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2397.29589647753        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      911
 NPARAMETR:  1.0408E+00  1.7726E+00  1.4393E+00  4.9152E-01  1.3683E+00  1.0287E+00  7.8464E-01  2.0474E+00  1.5797E+00  1.5245E+00
             2.0314E+00
 PARAMETER:  1.3998E-01  6.7246E-01  4.6418E-01 -6.1026E-01  4.1354E-01  1.2829E-01 -1.4253E-01  8.1659E-01  5.5724E-01  5.2169E-01
             8.0871E-01
 GRADIENT:  -7.9913E-01  1.1556E+01 -1.6251E+00  7.0686E+00  1.7666E+00 -3.9152E-01  8.1383E-02  1.1333E+00 -2.1366E+00 -5.3110E-02
            -3.7255E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2397.56248286728        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1089
 NPARAMETR:  1.0395E+00  1.9431E+00  9.3298E-01  3.7007E-01  1.3782E+00  1.0299E+00  7.4873E-01  1.3756E+00  1.9739E+00  1.5274E+00
             2.0326E+00
 PARAMETER:  1.3872E-01  7.6426E-01  3.0633E-02 -8.9406E-01  4.2080E-01  1.2943E-01 -1.8938E-01  4.1890E-01  7.8000E-01  5.2359E-01
             8.0930E-01
 GRADIENT:  -3.8714E+00  7.5207E+00 -7.6877E-01  6.0267E+00 -1.6597E-01 -1.4216E-01 -2.6453E-01  3.7405E-01  2.7171E+00 -1.4179E-01
             1.6177E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2397.67732795916        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1264
 NPARAMETR:  1.0418E+00  2.0835E+00  6.7560E-01  2.7784E-01  1.4129E+00  1.0304E+00  7.4970E-01  9.1596E-01  2.2465E+00  1.5604E+00
             2.0295E+00
 PARAMETER:  1.4098E-01  8.3407E-01 -2.9215E-01 -1.1807E+00  4.4563E-01  1.2995E-01 -1.8808E-01  1.2219E-02  9.0939E-01  5.4492E-01
             8.0777E-01
 GRADIENT:   5.2896E-01  1.2299E+01 -1.9439E-01  5.1176E+00 -2.0903E+00 -7.2964E-02  1.8116E-01  2.3781E-01  2.2875E+00 -3.5530E-01
            -3.7979E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2397.72919493372        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1439
 NPARAMETR:  1.0431E+00  2.1782E+00  5.3331E-01  2.1555E-01  1.4476E+00  1.0305E+00  7.4785E-01  6.1307E-01  2.4744E+00  1.5970E+00
             2.0273E+00
 PARAMETER:  1.4218E-01  8.7850E-01 -5.2864E-01 -1.4346E+00  4.6988E-01  1.3002E-01 -1.9055E-01 -3.8928E-01  1.0060E+00  5.6810E-01
             8.0670E-01
 GRADIENT:   2.8823E+00  1.3981E+01 -3.6862E-01  3.6576E+00 -1.8775E+00 -1.0305E-01  5.4712E-01  2.2225E-01  7.6153E-01 -4.9953E-01
            -1.5111E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2397.73737821851        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1614
 NPARAMETR:  1.0430E+00  2.2097E+00  4.9187E-01  1.9278E-01  1.4634E+00  1.0305E+00  7.4400E-01  5.0925E-01  2.5860E+00  1.6145E+00
             2.0269E+00
 PARAMETER:  1.4209E-01  8.9287E-01 -6.0953E-01 -1.5462E+00  4.8079E-01  1.3003E-01 -1.9572E-01 -5.7481E-01  1.0501E+00  5.7902E-01
             8.0648E-01
 GRADIENT:   2.7164E+00  1.0408E+01 -3.3509E-01  2.4473E+00 -1.1155E+00 -8.4940E-02  5.9866E-02  1.7896E-01  2.0469E-01 -2.9447E-01
            -1.4720E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2397.87214385751        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1795
 NPARAMETR:  1.0417E+00  2.2020E+00  4.9015E-01  1.8867E-01  1.4661E+00  1.0307E+00  7.4367E-01  1.9771E-01  2.5964E+00  1.6170E+00
             2.0281E+00
 PARAMETER:  1.4082E-01  8.8935E-01 -6.1305E-01 -1.5678E+00  4.8258E-01  1.3021E-01 -1.9616E-01 -1.5209E+00  1.0541E+00  5.8057E-01
             8.0710E-01
 GRADIENT:   3.3350E-01 -7.8257E+00  2.7316E-01 -5.0028E-01 -2.0086E-01  7.6956E-02  9.3478E-02  2.3843E-02 -5.7961E-01 -2.1208E-02
             4.0055E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2397.89111578570        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1979
 NPARAMETR:  1.0420E+00  2.1989E+00  4.8816E-01  1.8960E-01  1.4647E+00  1.0307E+00  7.4375E-01  2.7478E-02  2.6004E+00  1.6157E+00
             2.0276E+00
 PARAMETER:  1.4118E-01  8.8795E-01 -6.1711E-01 -1.5628E+00  4.8166E-01  1.3025E-01 -1.9606E-01 -3.4944E+00  1.0557E+00  5.7977E-01
             8.0685E-01
 GRADIENT:   1.1273E+00 -1.0665E+01  7.4911E-02 -4.1456E-01 -1.1427E-01  9.6251E-02  3.7066E-02  4.9147E-04  5.4875E-02  1.0053E-01
             6.6534E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2397.89427742392        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:     2079
 NPARAMETR:  1.0413E+00  2.1945E+00  4.8788E-01  1.9145E-01  1.4634E+00  1.0301E+00  7.4358E-01  1.0000E-02  2.5951E+00  1.6140E+00
             2.0272E+00
 PARAMETER:  1.4046E-01  8.8596E-01 -6.1768E-01 -1.5531E+00  4.8073E-01  1.2966E-01 -1.9628E-01 -5.1930E+00  1.0536E+00  5.7874E-01
             8.0663E-01
 GRADIENT:   1.8136E+02  5.5425E+02  1.3195E-01  2.1070E+01  2.0785E+01  1.8312E+01  6.7081E+00  0.0000E+00  1.0549E+01  7.2041E+00
             1.0992E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2397.89732468987        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2254
 NPARAMETR:  1.0421E+00  2.1841E+00  4.8801E-01  1.9502E-01  1.4635E+00  1.0308E+00  7.4457E-01  1.0000E-02  2.5719E+00  1.6132E+00
             2.0289E+00
 PARAMETER:  1.4115E-01  8.8695E-01 -6.1130E-01 -1.5498E+00  4.8047E-01  1.3021E-01 -1.9580E-01 -5.0753E+00  1.0517E+00  5.7845E-01
             8.0677E-01
 GRADIENT:  -1.5842E-02  2.1949E+00  2.1862E-02 -1.6907E-01 -3.0452E-02 -7.5638E-03 -2.7222E-02  0.0000E+00  1.1578E-01  7.4248E-03
            -1.3356E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2254
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.4462E-04 -3.2008E-02 -5.8799E-05  3.8336E-02 -2.8364E-02
 SE:             2.9628E-02  2.4915E-02  2.9965E-05  1.8135E-02  2.5390E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8533E-01  1.9891E-01  4.9730E-02  3.4523E-02  2.6394E-01

 ETASHRINKSD(%)  7.4218E-01  1.6531E+01  9.9900E+01  3.9245E+01  1.4939E+01
 ETASHRINKVR(%)  1.4788E+00  3.0329E+01  1.0000E+02  6.3089E+01  2.7646E+01
 EBVSHRINKSD(%)  9.8460E-01  1.3202E+01  9.9903E+01  5.0391E+01  1.0379E+01
 EBVSHRINKVR(%)  1.9595E+00  2.4661E+01  1.0000E+02  7.5390E+01  1.9680E+01
 RELATIVEINF(%)  9.7999E+01  1.4751E+01  5.3965E-05  4.4664E+00  5.1247E+01
 EPSSHRINKSD(%)  2.1984E+01
 EPSSHRINKVR(%)  3.9134E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          695
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1277.3245611544951     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2397.8973246898686     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1120.5727635353735     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.15
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2397.897       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  2.20E+00  4.91E-01  1.92E-01  1.46E+00  1.03E+00  7.44E-01  1.00E-02  2.59E+00  1.61E+00  2.03E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.01
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,      207.839
Stop Time:
Sat Oct 23 19:27:52 CDT 2021
