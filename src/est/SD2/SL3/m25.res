Sat Oct 23 19:14:20 CDT 2021
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
$DATA ../../../../data/SD2/SL3/dat25.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      799
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

 TOT. NO. OF OBS RECS:      699
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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2152.04416744249        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3816E+02 -1.2157E+02  2.4362E+01  9.9380E+01  7.7233E+01  5.6993E+01 -4.6021E+01 -6.3555E+01 -3.8597E+01  2.4929E+01
            -1.5093E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2539.33991508320        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0446E+00  1.3118E+00  9.7954E-01  8.7741E-01  1.0836E+00  9.9962E-01  1.0826E+00  9.8820E-01  1.0706E+00  8.3570E-01
             1.7307E+00
 PARAMETER:  1.4363E-01  3.7140E-01  7.9327E-02 -3.0775E-02  1.8031E-01  9.9620E-02  1.7933E-01  8.8134E-02  1.6823E-01 -7.9486E-02
             6.4851E-01
 GRADIENT:   3.0449E+02  1.0448E+02  1.4129E+01  3.6749E+01 -2.1564E+01  3.1487E+01  1.2573E+01 -1.7280E+00 -6.3118E+00 -1.0176E+01
            -1.8505E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2543.83886371730        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      245
 NPARAMETR:  1.0537E+00  1.4581E+00  9.3353E-01  8.1736E-01  1.1686E+00  9.2911E-01  8.9322E-01  1.0046E+00  1.1861E+00  1.0768E+00
             1.7458E+00
 PARAMETER:  1.5228E-01  4.7714E-01  3.1217E-02 -1.0167E-01  2.5583E-01  2.6473E-02 -1.2918E-02  1.0456E-01  2.7065E-01  1.7399E-01
             6.5723E-01
 GRADIENT:   1.3710E+02  1.9783E+01 -3.8556E+00  5.8382E+01 -2.0840E+01 -7.6230E+00 -3.9415E+00  1.5256E+00 -3.6190E+00  1.0681E+01
            -5.2252E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2557.82771440740        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      421
 NPARAMETR:  9.9504E-01  1.6204E+00  2.1845E+00  7.1308E-01  1.6791E+00  9.1284E-01  8.5545E-01  2.2239E+00  1.2993E+00  1.4584E+00
             1.7017E+00
 PARAMETER:  9.5025E-02  5.8265E-01  8.8138E-01 -2.3816E-01  6.1824E-01  8.8045E-03 -5.6130E-02  8.9926E-01  3.6185E-01  4.7733E-01
             6.3162E-01
 GRADIENT:  -4.1046E+00 -3.6304E+00 -1.0916E+01  3.0178E+01  3.5539E+01 -7.5983E+00  2.6984E+00 -2.2249E+00  3.3687E+00  6.0824E+00
            -1.7966E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2563.31547015306        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      597
 NPARAMETR:  9.9799E-01  1.9730E+00  2.1282E+00  4.6598E-01  1.7439E+00  9.2930E-01  7.7896E-01  3.2412E+00  1.5567E+00  1.5192E+00
             1.7102E+00
 PARAMETER:  9.7986E-02  7.7956E-01  8.5526E-01 -6.6362E-01  6.5614E-01  2.6681E-02 -1.4979E-01  1.2759E+00  5.4260E-01  5.1817E-01
             6.3663E-01
 GRADIENT:   2.3423E+00  9.5657E+00  3.0650E+00  4.1754E+00 -5.7820E+00 -2.7615E-01 -1.4557E+00 -2.0672E+00 -2.7195E+00  2.6343E+00
             1.5024E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2563.63998659333        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      777
 NPARAMETR:  9.9823E-01  2.0207E+00  2.0872E+00  4.2215E-01  1.7702E+00  9.3238E-01  7.7902E-01  3.5832E+00  1.6707E+00  1.5137E+00
             1.7084E+00
 PARAMETER:  9.8227E-02  8.0345E-01  8.3582E-01 -7.6239E-01  6.7108E-01  2.9985E-02 -1.4972E-01  1.3763E+00  6.1325E-01  5.1459E-01
             6.3556E-01
 GRADIENT:   3.1342E+00 -1.3475E+01  2.0259E+00 -1.0202E+00 -2.7829E+00  1.0463E+00  1.1890E+00  5.2471E-01  1.1286E+00 -1.1843E+00
             1.5878E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2563.84632275537        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      958
 NPARAMETR:  9.9538E-01  2.0359E+00  1.8395E+00  4.1124E-01  1.7841E+00  9.2972E-01  7.7509E-01  3.4568E+00  1.6502E+00  1.5230E+00
             1.7064E+00
 PARAMETER:  9.5372E-02  8.1093E-01  7.0948E-01 -7.8857E-01  6.7894E-01  2.7129E-02 -1.5478E-01  1.3403E+00  6.0087E-01  5.2068E-01
             6.3438E-01
 GRADIENT:  -4.1468E+00 -1.6142E+01 -3.5775E-01  4.0420E-01  5.8420E+00 -6.9874E-02 -3.7718E-01  1.7856E+00 -1.4838E+00 -3.5444E-01
            -5.7483E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2565.72507392685        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1139
 NPARAMETR:  1.0022E+00  2.2894E+00  4.9661E-01  2.7080E-01  1.6926E+00  9.2911E-01  7.4587E-01  1.4664E+00  2.0560E+00  1.5082E+00
             1.7095E+00
 PARAMETER:  1.0222E-01  9.2827E-01 -5.9996E-01 -1.2064E+00  6.2628E-01  2.6473E-02 -1.9321E-01  4.8279E-01  8.2075E-01  5.1095E-01
             6.3623E-01
 GRADIENT:   1.0838E+01  6.3256E+01 -2.7830E+00  1.6647E+01 -3.2701E+00 -9.0597E-01  3.6752E-01  2.4378E+00 -7.2510E+00  9.1654E+00
            -1.9932E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2569.05669832448        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1315
 NPARAMETR:  9.9802E-01  2.3879E+00  3.1030E-01  1.5912E-01  1.7504E+00  9.3151E-01  7.2312E-01  8.4069E-01  2.7784E+00  1.4789E+00
             1.6992E+00
 PARAMETER:  9.8014E-02  9.7042E-01 -1.0702E+00 -1.7381E+00  6.5983E-01  2.9048E-02 -2.2418E-01 -7.3538E-02  1.1219E+00  4.9131E-01
             6.3014E-01
 GRADIENT:   1.3534E+00 -1.9289E+01 -6.2637E-01 -1.5676E+00 -6.8229E+00  3.4169E-01  7.6781E-01  1.2619E+00 -1.9563E-01 -4.0935E+00
            -4.0754E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2569.09895421102        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1490
 NPARAMETR:  9.9669E-01  2.4200E+00  2.7731E-01  1.4223E-01  1.7791E+00  9.3093E-01  7.2014E-01  7.6783E-01  2.9316E+00  1.5111E+00
             1.6987E+00
 PARAMETER:  9.6687E-02  9.8375E-01 -1.1826E+00 -1.8503E+00  6.7609E-01  2.8433E-02 -2.2830E-01 -1.6419E-01  1.1756E+00  5.1284E-01
             6.2988E-01
 GRADIENT:  -2.1246E+00 -1.2020E+01 -1.6222E+00 -1.6176E-02 -2.1284E+00  9.7090E-02  5.9723E-01  1.0801E+00  7.7831E-01 -2.0672E+00
            -3.8812E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2569.10158598136        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1665
 NPARAMETR:  9.9699E-01  2.4299E+00  2.7511E-01  1.3751E-01  1.7895E+00  9.3083E-01  7.1858E-01  7.5413E-01  2.9808E+00  1.5236E+00
             1.7002E+00
 PARAMETER:  9.6987E-02  9.8784E-01 -1.1906E+00 -1.8840E+00  6.8195E-01  2.8326E-02 -2.3047E-01 -1.8220E-01  1.1922E+00  5.2109E-01
             6.3074E-01
 GRADIENT:  -1.4517E+00 -8.0613E+00 -1.3248E+00  3.8508E-02 -1.0188E+00  7.0446E-02  2.5651E-01  1.0294E+00  8.9803E-01 -1.2412E+00
            -2.5494E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2569.10498485784        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1840
 NPARAMETR:  9.9729E-01  2.4385E+00  2.7606E-01  1.3356E-01  1.7993E+00  9.3075E-01  7.1710E-01  7.3140E-01  3.0235E+00  1.5354E+00
             1.7018E+00
 PARAMETER:  9.7290E-02  9.9138E-01 -1.1871E+00 -1.9132E+00  6.8742E-01  2.8230E-02 -2.3254E-01 -2.1280E-01  1.2064E+00  5.2882E-01
             6.3171E-01
 GRADIENT:  -7.6474E-01 -4.1503E+00 -9.2346E-01 -3.4239E-03 -1.4126E-01  4.9170E-02 -6.6732E-02  9.6229E-01  9.5288E-01 -4.9686E-01
            -1.2029E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2569.11717606508        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2016
 NPARAMETR:  9.9778E-01  2.4497E+00  2.8260E-01  1.2876E-01  1.8139E+00  9.3061E-01  7.1496E-01  6.6466E-01  3.0742E+00  1.5531E+00
             1.7046E+00
 PARAMETER:  9.7776E-02  9.9596E-01 -1.1637E+00 -1.9498E+00  6.9546E-01  2.8083E-02 -2.3552E-01 -3.0848E-01  1.2230E+00  5.4026E-01
             6.3335E-01
 GRADIENT:   3.4443E-01  1.5264E+00 -2.4687E-01 -1.1206E-01  1.0078E+00  1.9405E-02 -5.2027E-01  8.0027E-01  9.5512E-01  5.4089E-01
             7.6256E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2571.10525794226        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2201
 NPARAMETR:  9.9668E-01  2.5007E+00  3.3453E-01  1.0554E-01  1.8789E+00  9.2941E-01  7.1790E-01  2.1460E-01  3.1711E+00  1.6241E+00
             1.7075E+00
 PARAMETER:  9.6675E-02  1.0166E+00 -9.9504E-01 -2.1487E+00  7.3068E-01  2.6798E-02 -2.3142E-01 -1.4390E+00  1.2541E+00  5.8497E-01
             6.3503E-01
 GRADIENT:  -2.5722E+00  3.2496E+01  1.7186E+00 -6.2033E-01  7.5234E+00 -7.8621E-01  3.0206E+00  3.5844E-02 -1.4009E+01  2.3739E+00
             6.4308E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2573.21090377168        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2377
 NPARAMETR:  9.9777E-01  2.5654E+00  2.5612E-01  5.4606E-02  1.9165E+00  9.3035E-01  7.1446E-01  8.2764E-02  4.7119E+00  1.6651E+00
             1.6980E+00
 PARAMETER:  9.7767E-02  1.0421E+00 -1.2621E+00 -2.8076E+00  7.5051E-01  2.7810E-02 -2.3623E-01 -2.3918E+00  1.6501E+00  6.0987E-01
             6.2947E-01
 GRADIENT:   2.3304E-01  2.3719E+01  1.4271E+00 -3.2110E+00  9.1394E+00 -2.0158E-01  3.8918E+00  1.1276E-02 -7.8066E+00  6.5524E+00
             5.5404E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2573.39162325141        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     2513
 NPARAMETR:  9.9718E-01  2.5658E+00  2.4563E-01  4.9330E-02  1.9012E+00  9.3013E-01  7.1239E-01  5.8514E-02  4.9230E+00  1.6324E+00
             1.6954E+00
 PARAMETER:  9.7180E-02  1.0423E+00 -1.3039E+00 -2.9092E+00  7.4251E-01  2.7564E-02 -2.3913E-01 -2.7385E+00  1.6939E+00  5.9005E-01
             6.2790E-01
 GRADIENT:  -9.8682E-01  1.2155E+01  5.9100E-01 -2.9688E+00  4.2541E+00 -1.7044E-01  2.5898E+00  5.9314E-03 -6.0566E+00  1.8584E+00
             2.2117E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2573.42108385577        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2705
 NPARAMETR:  9.9766E-01  2.5660E+00  2.4515E-01  4.9546E-02  1.8954E+00  9.3030E-01  7.1537E-01  2.1696E-02  4.9370E+00  1.6266E+00
             1.6940E+00
 PARAMETER:  9.7654E-02  1.0423E+00 -1.3059E+00 -2.9044E+00  7.3947E-01  2.7761E-02 -2.3495E-01 -3.6937E+00  1.6970E+00  5.8649E-01
             6.2710E-01
 GRADIENT:  -9.3627E-02 -1.5775E+01  2.2864E-02  1.1577E+01  6.2595E-01  5.2248E-02  2.5796E-01  8.6143E-04  2.2996E+01 -7.2918E-02
            -6.1321E-02
 NUMSIGDIG:         3.4         2.9         3.3         2.4         3.0         2.9         3.0         0.6         2.4         3.4
                    4.2

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2705
 NO. OF SIG. DIGITS IN FINAL EST.:  0.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.6144E-04 -2.8363E-02  1.6569E-06  3.5479E-02 -2.8573E-02
 SE:             2.9674E-02  2.7163E-02  2.5956E-05  1.4120E-02  2.5449E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9297E-01  2.9639E-01  9.4910E-01  1.1981E-02  2.6154E-01

 ETASHRINKSD(%)  5.8943E-01  9.0021E+00  9.9913E+01  5.2696E+01  1.4744E+01
 ETASHRINKVR(%)  1.1754E+00  1.7194E+01  1.0000E+02  7.7624E+01  2.7313E+01
 EBVSHRINKSD(%)  8.8517E-01  5.1079E+00  9.9780E+01  6.8039E+01  9.6270E+00
 EBVSHRINKVR(%)  1.7625E+00  9.9549E+00  1.0000E+02  8.9785E+01  1.8327E+01
 RELATIVEINF(%)  9.8216E+01  5.3065E+01  4.4281E-04  5.8737E+00  7.6514E+01
 EPSSHRINKSD(%)  2.1920E+01
 EPSSHRINKVR(%)  3.9035E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          699
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1284.6760694201323     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2573.4210838557747     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1288.7450144356424     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.99
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2573.421       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.98E-01  2.57E+00  2.45E-01  4.96E-02  1.90E+00  9.30E-01  7.15E-01  2.25E-02  4.94E+00  1.63E+00  1.69E+00
 


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
 #CPUT: Total CPU Time in Seconds,      250.833
Stop Time:
Sat Oct 23 19:14:55 CDT 2021
