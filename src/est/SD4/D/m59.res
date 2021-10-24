Sun Oct 24 04:21:27 CDT 2021
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
$DATA ../../../../data/SD4/D/dat59.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m59.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1679.55626626739        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6683E+02  5.7213E+00 -2.2869E+01  4.2326E+01  1.3448E+01  1.6597E+01 -1.3330E+01  4.7884E+00 -2.0704E+01  1.9270E+01
             1.8393E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1682.39363593475        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      191
 NPARAMETR:  1.1240E+00  1.0697E+00  1.1579E+00  1.0338E+00  1.0642E+00  1.3909E+00  1.1817E+00  9.9244E-01  1.2344E+00  8.6390E-01
             9.3581E-01
 PARAMETER:  2.1685E-01  1.6738E-01  2.4657E-01  1.3325E-01  1.6226E-01  4.2997E-01  2.6693E-01  9.2414E-02  3.1056E-01 -4.6294E-02
             3.3661E-02
 GRADIENT:   7.9610E+01  2.8679E+01  1.7163E+01  4.0022E+01 -8.9422E+00  6.2467E+01  1.1262E+00 -9.9980E+00  2.1477E+01 -1.0455E+01
            -1.5405E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1685.60169843005        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  1.1058E+00  8.7564E-01  1.4166E+00  1.1777E+00  1.0540E+00  1.3353E+00  1.4034E+00  1.5258E+00  1.0892E+00  8.5148E-01
             9.3017E-01
 PARAMETER:  2.0059E-01 -3.2798E-02  4.4827E-01  2.6353E-01  1.5261E-01  3.8915E-01  4.3887E-01  5.2250E-01  1.8540E-01 -6.0784E-02
             2.7613E-02
 GRADIENT:   6.6388E+01  3.3481E+01  5.7103E+00  5.3121E+01 -2.9915E+01  5.2607E+01  6.7994E+00  8.3498E+00  1.3295E+01 -5.3014E+00
            -1.4472E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1694.47555800267        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      544
 NPARAMETR:  1.0536E+00  8.8034E-01  1.4069E+00  1.1249E+00  1.0899E+00  1.1379E+00  1.2173E+00  1.3522E+00  1.0608E+00  9.5015E-01
             9.6173E-01
 PARAMETER:  1.5219E-01 -2.7442E-02  4.4140E-01  2.1766E-01  1.8611E-01  2.2917E-01  2.9668E-01  4.0171E-01  1.5901E-01  4.8860E-02
             6.0974E-02
 GRADIENT:   5.8160E-01  3.0077E+00  7.1044E-01 -1.0979E+00 -3.0391E+00 -1.8577E+00  8.0350E-01  1.4748E-02  1.8652E-01  6.4784E-01
             2.1727E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1694.83619049356        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      720
 NPARAMETR:  1.0485E+00  6.3405E-01  1.7011E+00  1.3001E+00  1.0857E+00  1.1498E+00  1.2151E+00  1.4934E+00  9.9651E-01  9.6439E-01
             9.6336E-01
 PARAMETER:  1.4732E-01 -3.5563E-01  6.3126E-01  3.6246E-01  1.8221E-01  2.3956E-01  2.9481E-01  5.0106E-01  9.6501E-02  6.3737E-02
             6.2671E-02
 GRADIENT:  -3.1272E+00  9.3296E+00  7.0757E+00  1.4900E+01 -8.0211E+00  3.5144E+00  1.1650E-01 -2.9270E+00  1.1984E+00 -9.9688E-01
             4.1458E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1695.14654406833        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      898
 NPARAMETR:  1.0487E+00  4.3350E-01  1.8034E+00  1.4291E+00  1.0567E+00  1.1358E+00  1.0318E+00  1.5863E+00  9.4095E-01  9.6462E-01
             9.6031E-01
 PARAMETER:  1.4759E-01 -7.3587E-01  6.8966E-01  4.5706E-01  1.5514E-01  2.2736E-01  1.3130E-01  5.6138E-01  3.9131E-02  6.3980E-02
             5.9498E-02
 GRADIENT:   1.3522E+00  6.7716E+00  1.2087E+00  1.6719E+01 -5.8390E+00 -6.1249E-01 -5.1112E-01  4.7830E-01 -1.1244E+00  3.0267E-01
            -8.2830E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1695.49563700720        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  1.0445E+00  2.9088E-01  1.8853E+00  1.5057E+00  1.0423E+00  1.1312E+00  1.1072E+00  1.6508E+00  9.0222E-01  9.5747E-01
             9.6102E-01
 PARAMETER:  1.4350E-01 -1.1348E+00  7.3408E-01  5.0929E-01  1.4140E-01  2.2329E-01  2.0188E-01  6.0126E-01 -2.8958E-03  5.6538E-02
             6.0240E-02
 GRADIENT:   7.6187E+02  5.5163E+01  9.7980E+00  1.0408E+03  9.5749E+00  1.7648E+02  7.8165E-01  4.5095E+00  2.0449E+01  1.0463E+00
             6.8982E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1695.57449968505        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1220             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0481E+00  2.7335E-01  1.8897E+00  1.5121E+00  1.0393E+00  1.1391E+00  1.4244E+00  1.6339E+00  8.9010E-01  9.5552E-01
             9.6188E-01
 PARAMETER:  1.4696E-01 -1.1970E+00  7.3644E-01  5.1351E-01  1.3856E-01  2.3021E-01  4.5377E-01  5.9100E-01 -1.6421E-02  5.4501E-02
             6.1132E-02
 GRADIENT:   7.8499E+02  5.1049E+01  1.0887E+01  1.0602E+03  1.1277E+01  1.8255E+02  2.7313E+00  3.1536E+00  2.0019E+01  9.8586E-01
             9.8857E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1695.59743619866        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1400             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0512E+00  2.6348E-01  1.8909E+00  1.5198E+00  1.0342E+00  1.1429E+00  1.4870E+00  1.6344E+00  8.8227E-01  9.5104E-01
             9.6172E-01
 PARAMETER:  1.4995E-01 -1.2338E+00  7.3703E-01  5.1859E-01  1.3364E-01  2.3354E-01  4.9676E-01  5.9130E-01 -2.5254E-02  4.9800E-02
             6.0965E-02
 GRADIENT:   8.0616E+02  4.9749E+01  1.1952E+01  1.0840E+03  8.1882E+00  1.8498E+02  2.9440E+00  3.2111E+00  1.8569E+01  1.0091E+00
             8.0073E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1695.61038146480        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1559
 NPARAMETR:  1.0486E+00  2.5875E-01  1.8899E+00  1.5274E+00  1.0334E+00  1.1390E+00  1.4662E+00  1.6339E+00  8.7986E-01  9.4936E-01
             9.6176E-01
 PARAMETER:  1.4742E-01 -1.2519E+00  7.3651E-01  5.2360E-01  1.3284E-01  2.3014E-01  4.8267E-01  5.9096E-01 -2.7988E-02  4.8037E-02
             6.1009E-02
 GRADIENT:   5.5079E+00  1.2720E+00  1.6779E-01 -1.3414E+01  4.7593E-01  1.2690E+00 -1.1072E-02  1.9452E-02  2.8210E-01  1.2439E-01
            -9.1301E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1695.64037212172        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1742             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0477E+00  2.3271E-01  1.8929E+00  1.5361E+00  1.0298E+00  1.1385E+00  1.5565E+00  1.6352E+00  8.7438E-01  9.4674E-01
             9.6200E-01
 PARAMETER:  1.4660E-01 -1.3580E+00  7.3810E-01  5.2922E-01  1.2936E-01  2.2970E-01  5.4245E-01  5.9176E-01 -3.4243E-02  4.5268E-02
             6.1255E-02
 GRADIENT:   7.8188E+02  4.2929E+01  1.0470E+01  1.1299E+03  1.2087E+01  1.8120E+02  2.7764E+00  3.1402E+00  1.9130E+01  5.8807E-01
             9.4175E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1695.67914935825        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1924             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0471E+00  2.1589E-01  1.9004E+00  1.5492E+00  1.0257E+00  1.1382E+00  1.5665E+00  1.6379E+00  8.6886E-01  9.4666E-01
             9.6199E-01
 PARAMETER:  1.4601E-01 -1.4330E+00  7.4207E-01  5.3774E-01  1.2535E-01  2.2944E-01  5.4882E-01  5.9339E-01 -4.0572E-02  4.5179E-02
             6.1253E-02
 GRADIENT:   7.7771E+02  3.9873E+01  1.1191E+01  1.1693E+03  9.7415E+00  1.8088E+02  2.4371E+00  3.0405E+00  1.8987E+01  9.2146E-01
             9.1449E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1695.69038638900        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2108
 NPARAMETR:  1.0458E+00  2.1043E-01  1.9011E+00  1.5549E+00  1.0231E+00  1.1377E+00  1.5613E+00  1.6402E+00  8.6531E-01  9.4263E-01
             9.6185E-01
 PARAMETER:  1.4479E-01 -1.4586E+00  7.4245E-01  5.4140E-01  1.2280E-01  2.2899E-01  5.4550E-01  5.9479E-01 -4.4668E-02  4.0922E-02
             6.1102E-02
 GRADIENT:   2.1520E+00  7.3740E-01  9.8668E-02 -2.0673E+01  6.0919E-01  1.0217E+00 -1.3210E-02 -2.8631E-02  8.4898E-02  5.6954E-02
            -4.5165E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1695.72475650960        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2291             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0466E+00  1.8343E-01  1.9061E+00  1.5687E+00  1.0189E+00  1.1378E+00  1.6586E+00  1.6422E+00  8.5982E-01  9.4113E-01
             9.6197E-01
 PARAMETER:  1.4552E-01 -1.5959E+00  7.4503E-01  5.5022E-01  1.1869E-01  2.2912E-01  6.0596E-01  5.9606E-01 -5.1027E-02  3.9324E-02
             6.1226E-02
 GRADIENT:   7.7411E+02  3.2855E+01  1.0959E+01  1.2286E+03  1.0044E+01  1.8007E+02  2.1897E+00  3.1013E+00  1.9343E+01  7.3317E-01
             8.7989E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1695.73818670575        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2475
 NPARAMETR:  1.0463E+00  1.7896E-01  1.9085E+00  1.5738E+00  1.0162E+00  1.1376E+00  1.6668E+00  1.6443E+00  8.5586E-01  9.3877E-01
             9.6190E-01
 PARAMETER:  1.4522E-01 -1.6206E+00  7.4634E-01  5.5347E-01  1.1603E-01  2.2890E-01  6.1089E-01  5.9729E-01 -5.5648E-02  3.6810E-02
             6.1153E-02
 GRADIENT:   3.7835E+00  5.9237E-01  1.7352E-01 -2.3895E+01  1.7323E-01  1.1195E+00 -7.1420E-03 -5.1644E-02  2.6420E-02  1.3378E-01
            -3.5188E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1695.76732558958        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2660
 NPARAMETR:  1.0460E+00  1.5405E-01  1.9143E+00  1.5878E+00  1.0126E+00  1.1374E+00  1.7792E+00  1.6472E+00  8.5081E-01  9.3740E-01
             9.6203E-01
 PARAMETER:  1.4501E-01 -1.7705E+00  7.4934E-01  5.6232E-01  1.1247E-01  2.2874E-01  6.7615E-01  5.9908E-01 -6.1564E-02  3.5355E-02
             6.1292E-02
 GRADIENT:   4.1107E+00  2.4862E-01 -4.9268E-01 -2.8463E+01  1.8949E+00  1.1624E+00  1.9445E-02 -1.0947E-01  9.7569E-01  1.5415E-01
             7.3446E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1695.77213123088        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2840
 NPARAMETR:  1.0436E+00  1.3570E-01  1.9195E+00  1.6048E+00  1.0094E+00  1.1342E+00  1.8047E+00  1.6502E+00  8.4593E-01  9.3578E-01
             9.6203E-01
 PARAMETER:  1.4264E-01 -1.8973E+00  7.5206E-01  5.7300E-01  1.0938E-01  2.2589E-01  6.9038E-01  6.0090E-01 -6.7319E-02  3.3622E-02
             6.1295E-02
 GRADIENT:   3.4278E-01  7.7824E-01 -9.6013E-01 -1.9640E+01  1.7748E+00  7.1413E-02 -1.8133E-03 -2.0830E-01  9.7159E-01  1.3047E-01
            -1.3654E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1695.80562543629        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     3023
 NPARAMETR:  1.0447E+00  1.2175E-01  1.9261E+00  1.6107E+00  1.0058E+00  1.1359E+00  1.8630E+00  1.6547E+00  8.4015E-01  9.3321E-01
             9.6203E-01
 PARAMETER:  1.4369E-01 -2.0058E+00  7.5550E-01  5.7666E-01  1.0574E-01  2.2744E-01  7.2217E-01  6.0363E-01 -7.4180E-02  3.0876E-02
             6.1293E-02
 GRADIENT:   2.6422E+00  5.3968E-01 -1.3060E-01 -2.5674E+01  2.0466E-01  7.7632E-01 -7.8347E-03 -1.6585E-01  1.8495E-01  2.0158E-01
            -2.4283E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1695.81458653233        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     3206
 NPARAMETR:  1.0443E+00  1.1134E-01  1.9286E+00  1.6170E+00  1.0049E+00  1.1354E+00  2.0646E+00  1.6574E+00  8.3856E-01  9.3162E-01
             9.6210E-01
 PARAMETER:  1.4330E-01 -2.0951E+00  7.5677E-01  5.8055E-01  1.0493E-01  2.2697E-01  8.2492E-01  6.0523E-01 -7.6069E-02  2.9172E-02
             6.1362E-02
 GRADIENT:   2.2613E+00  4.1176E-01 -7.8923E-01 -2.7174E+01  1.8374E+00  6.2588E-01  1.8282E-02 -1.4850E-01  1.0613E+00  4.4807E-02
             2.0831E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1695.82275956449        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     3395
 NPARAMETR:  1.0465E+00  9.3049E-02  1.9394E+00  1.6225E+00  1.0008E+00  1.1387E+00  2.1905E+00  1.6635E+00  8.3201E-01  9.3119E-01
             9.6218E-01
 PARAMETER:  1.4546E-01 -2.2746E+00  7.6240E-01  5.8399E-01  1.0082E-01  2.2992E-01  8.8412E-01  6.0890E-01 -8.3907E-02  2.8706E-02
             6.1445E-02
 GRADIENT:   6.6425E+00  2.1223E-02  4.9418E-01 -3.9364E+01 -2.3082E-01  1.9078E+00  1.6095E-02 -1.0557E-01  3.9512E-01  4.3519E-01
             1.2715E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1695.84322949898        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3577
 NPARAMETR:  1.0448E+00  8.6460E-02  1.9402E+00  1.6305E+00  1.0005E+00  1.1365E+00  2.3048E+00  1.6665E+00  8.3021E-01  9.2726E-01
             9.6211E-01
 PARAMETER:  1.4386E-01 -2.3481E+00  7.6281E-01  5.8887E-01  1.0055E-01  2.2793E-01  9.3500E-01  6.1073E-01 -8.6081E-02  2.4473E-02
             6.1374E-02
 GRADIENT:   3.9786E+00  2.3142E-01 -4.4786E-01 -3.3024E+01  1.3268E+00  1.1296E+00  1.4250E-02 -1.0440E-01  4.8822E-01 -7.7393E-02
            -2.0839E-02

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1695.85147924857        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3762
 NPARAMETR:  1.0455E+00  7.0289E-02  1.9498E+00  1.6383E+00  9.9710E-01  1.1376E+00  2.4011E+00  1.6715E+00  8.2488E-01  9.2801E-01
             9.6224E-01
 PARAMETER:  1.4449E-01 -2.5551E+00  7.6775E-01  5.9367E-01  9.7092E-02  2.2893E-01  9.7594E-01  6.1375E-01 -9.2519E-02  2.5289E-02
             6.1506E-02
 GRADIENT:   5.5480E+00  1.1653E-01  5.0769E-01 -3.8524E+01 -7.3231E-01  1.6088E+00  7.3653E-03 -8.4253E-02  2.7007E-02  3.6816E-01
             8.4762E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1695.86528060272        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3944             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0443E+00  5.8582E-02  1.9526E+00  1.6479E+00  9.9655E-01  1.1361E+00  2.9169E+00  1.6744E+00  8.2304E-01  9.2587E-01
             9.6224E-01
 PARAMETER:  1.4337E-01 -2.7373E+00  7.6914E-01  5.9953E-01  9.6546E-02  2.2757E-01  1.1705E+00  6.1544E-01 -9.4749E-02  2.2981E-02
             6.1514E-02
 GRADIENT:   7.5863E+02  8.2296E+00  1.1918E+01  1.4879E+03  7.6170E+00  1.7676E+02  1.3234E+00  3.3340E+00  1.9928E+01  7.0025E-01
             8.4619E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1695.86939683389        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     4123
 NPARAMETR:  1.0435E+00  5.5010E-02  1.9558E+00  1.6507E+00  9.9512E-01  1.1359E+00  2.6214E+00  1.6772E+00  8.2070E-01  9.2530E-01
             9.6223E-01
 PARAMETER:  1.4256E-01 -2.8002E+00  7.7082E-01  6.0121E-01  9.5109E-02  2.2744E-01  1.0637E+00  6.1710E-01 -9.7596E-02  2.2365E-02
             6.1499E-02
 GRADIENT:   2.5504E+00  2.0095E-01  3.9831E-02 -3.4678E+01 -2.9212E-01  1.0822E+00  4.2829E-03 -6.7606E-02  5.8664E-02  1.7899E-01
             1.5994E-02

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1695.87235295744        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     4305
 NPARAMETR:  1.0445E+00  5.2036E-02  1.9581E+00  1.6526E+00  9.9415E-01  1.1366E+00  2.9142E+00  1.6795E+00  8.1967E-01  9.2477E-01
             9.6228E-01
 PARAMETER:  1.4358E-01 -2.8558E+00  7.7199E-01  6.0238E-01  9.4129E-02  2.2801E-01  1.1696E+00  6.1850E-01 -9.8856E-02  2.1793E-02
             6.1554E-02
 GRADIENT:   4.4317E+00  2.1033E-01  3.5524E-01 -3.4932E+01 -1.2003E+00  1.3121E+00  1.1466E-02  1.3030E-02  1.2057E-01  2.7780E-01
             4.5860E-02

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1695.87392019754        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     4484
 NPARAMETR:  1.0435E+00  4.8858E-02  1.9574E+00  1.6542E+00  9.9450E-01  1.1356E+00  2.8413E+00  1.6797E+00  8.1945E-01  9.2321E-01
             9.6224E-01
 PARAMETER:  1.4262E-01 -2.9188E+00  7.7163E-01  6.0329E-01  9.4487E-02  2.2712E-01  1.1443E+00  6.1862E-01 -9.9125E-02  2.0099E-02
             6.1508E-02
 GRADIENT:   2.8473E+00  1.5502E-01 -3.1118E-01 -3.6066E+01  6.0338E-01  9.7103E-01  7.7386E-03 -2.5781E-02  3.0149E-01 -3.2159E-02
             9.7507E-05

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1695.87531666408        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     4646
 NPARAMETR:  1.0444E+00  4.8182E-02  1.9591E+00  1.6563E+00  9.9407E-01  1.1365E+00  2.9346E+00  1.6806E+00  8.1877E-01  9.2346E-01
             9.6226E-01
 PARAMETER:  1.4341E-01 -2.9328E+00  7.7248E-01  6.0458E-01  9.4049E-02  2.2794E-01  1.1766E+00  6.1913E-01 -9.9948E-02  2.0371E-02
             6.1530E-02
 GRADIENT:   4.2447E+00  2.2509E-01 -5.4793E-02 -3.3023E+01 -3.8561E-01  1.3049E+00  8.0618E-03 -3.0935E-02  1.4432E-01  7.6010E-02
            -1.7850E-02

0ITERATION NO.:  135    OBJECTIVE VALUE:  -1695.87911727175        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     4833
 NPARAMETR:  1.0445E+00  4.2002E-02  1.9614E+00  1.6584E+00  9.9404E-01  1.1365E+00  3.1498E+00  1.6823E+00  8.1769E-01  9.2305E-01
             9.6233E-01
 PARAMETER:  1.4354E-01 -3.0700E+00  7.7365E-01  6.0587E-01  9.4024E-02  2.2794E-01  1.2473E+00  6.2014E-01 -1.0128E-01  1.9931E-02
             6.1602E-02
 GRADIENT:   4.6651E+00  1.2844E-01 -4.5187E-01 -3.6973E+01  9.8495E-01  1.3350E+00  1.0904E-02 -4.9213E-02  4.3924E-01 -2.5051E-02
             2.7582E-02

0ITERATION NO.:  140    OBJECTIVE VALUE:  -1695.88069307167        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     5013
 NPARAMETR:  1.0437E+00  3.3411E-02  1.9661E+00  1.6672E+00  9.9325E-01  1.1355E+00  3.2094E+00  1.6845E+00  8.1602E-01  9.2336E-01
             9.6238E-01
 PARAMETER:  1.4277E-01 -3.2989E+00  7.7604E-01  6.1116E-01  9.3223E-02  2.2703E-01  1.2661E+00  6.2148E-01 -1.0332E-01  2.0268E-02
             6.1658E-02
 GRADIENT:   3.4951E+00  1.7503E-01 -6.5451E-01 -3.1440E+01  8.2492E-01  9.7939E-01  5.1070E-03 -1.4460E-01  6.8310E-01  3.5786E-02
             2.7479E-03

0ITERATION NO.:  145    OBJECTIVE VALUE:  -1695.88982011620        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     5203
 NPARAMETR:  1.0439E+00  2.6150E-02  1.9708E+00  1.6688E+00  9.9199E-01  1.1359E+00  4.1736E+00  1.6882E+00  8.1322E-01  9.2235E-01
             9.6235E-01
 PARAMETER:  1.4293E-01 -3.5439E+00  7.7842E-01  6.1211E-01  9.1960E-02  2.2742E-01  1.5288E+00  6.2364E-01 -1.0675E-01  1.9171E-02
             6.1622E-02
 GRADIENT:   4.0289E+00  1.0120E-01 -2.8964E-01 -3.7838E+01  5.5591E-01  1.1868E+00  1.2515E-02 -8.2874E-02  4.9195E-01  7.6991E-02
             2.4182E-02

0ITERATION NO.:  150    OBJECTIVE VALUE:  -1695.89147079157        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:     5320
 NPARAMETR:  1.0442E+00  2.3721E-02  1.9740E+00  1.6701E+00  9.9127E-01  1.1365E+00  4.4840E+00  1.6908E+00  8.1175E-01  9.2192E-01
             9.6237E-01
 PARAMETER:  1.4325E-01 -3.6414E+00  7.8008E-01  6.1289E-01  9.1230E-02  2.2795E-01  1.6005E+00  6.2521E-01 -1.0856E-01  1.8702E-02
             6.1648E-02
 GRADIENT:   7.5790E+02  2.8968E+00  1.2825E+01  1.5649E+03  5.8323E+00  1.7667E+02  7.0288E-01  3.5739E+00  2.1022E+01  7.6196E-01
             8.2973E-01

0ITERATION NO.:  155    OBJECTIVE VALUE:  -1695.89175674026        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     5498
 NPARAMETR:  1.0440E+00  2.1446E-02  1.9723E+00  1.6714E+00  9.9172E-01  1.1362E+00  3.9973E+00  1.6917E+00  8.1126E-01  9.1915E-01
             9.6230E-01
 PARAMETER:  1.4305E-01 -3.7422E+00  7.7921E-01  6.1364E-01  9.1682E-02  2.2765E-01  1.4856E+00  6.2573E-01 -1.0917E-01  1.5692E-02
             6.1573E-02
 GRADIENT:   4.3863E+00  6.5061E-02 -7.1409E-01 -3.9023E+01  1.5998E+00  1.3129E+00  4.8444E-03 -2.0910E-02  6.1282E-02 -3.0844E-01
            -6.5165E-02

0ITERATION NO.:  160    OBJECTIVE VALUE:  -1695.89420760060        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     5682             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0439E+00  1.8870E-02  1.9742E+00  1.6730E+00  9.9092E-01  1.1361E+00  4.3722E+00  1.6923E+00  8.1080E-01  9.1986E-01
             9.6232E-01
 PARAMETER:  1.4300E-01 -3.8702E+00  7.8015E-01  6.1464E-01  9.0878E-02  2.2762E-01  1.5753E+00  6.2608E-01 -1.0974E-01  1.6471E-02
             6.1590E-02
 GRADIENT:   7.5614E+02  1.9410E+00  1.2345E+01  1.5752E+03  7.0081E+00  1.7628E+02  4.3608E-01  3.6051E+00  2.1335E+01  4.9145E-01
             7.6913E-01

0ITERATION NO.:  165    OBJECTIVE VALUE:  -1695.89489983772        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     5839
 NPARAMETR:  1.0438E+00  1.8198E-02  1.9766E+00  1.6742E+00  9.9008E-01  1.1359E+00  4.8759E+00  1.6931E+00  8.1025E-01  9.2125E-01
             9.6238E-01
 PARAMETER:  1.4283E-01 -3.9064E+00  7.8136E-01  6.1536E-01  9.0034E-02  2.2747E-01  1.6843E+00  6.2653E-01 -1.1041E-01  1.7981E-02
             6.1653E-02
 GRADIENT:   7.5492E+02  2.0812E+00  1.2990E+01  1.5794E+03  5.2634E+00  1.7616E+02  5.2307E-01  3.6397E+00  2.1401E+01  8.0473E-01
             8.2476E-01

0ITERATION NO.:  170    OBJECTIVE VALUE:  -1695.89594168004        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     6026
 NPARAMETR:  1.0435E+00  1.5308E-02  1.9762E+00  1.6771E+00  9.9051E-01  1.1355E+00  4.5841E+00  1.6933E+00  8.1007E-01  9.2011E-01
             9.6237E-01
 PARAMETER:  1.4258E-01 -4.0794E+00  7.8115E-01  6.1707E-01  9.0465E-02  2.2711E-01  1.6226E+00  6.2669E-01 -1.1064E-01  1.6736E-02
             6.1640E-02
 GRADIENT:   3.7848E+00  7.0033E-02 -4.4293E-01 -3.6088E+01  4.9971E-01  1.1502E+00  3.9619E-03 -4.3788E-02  3.1439E-01 -6.3506E-02
            -2.9690E-02

0ITERATION NO.:  175    OBJECTIVE VALUE:  -1695.89810927744        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     6192
 NPARAMETR:  1.0433E+00  1.2027E-02  1.9781E+00  1.6783E+00  9.9016E-01  1.1354E+00  6.0792E+00  1.6943E+00  8.0927E-01  9.2030E-01
             9.6241E-01
 PARAMETER:  1.4242E-01 -4.3206E+00  7.8214E-01  6.1777E-01  9.0109E-02  2.2694E-01  1.9049E+00  6.2725E-01 -1.1162E-01  1.6950E-02
             6.1685E-02
 GRADIENT:   3.4981E+00  5.9756E-02 -3.7330E-01 -3.8108E+01  5.5079E-01  1.0568E+00  7.9060E-03 -5.1538E-02  4.4100E-01 -4.2917E-03
             1.2268E-02

0ITERATION NO.:  180    OBJECTIVE VALUE:  -1695.89950524053        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     6381
 NPARAMETR:  1.0438E+00  1.0000E-02  1.9808E+00  1.6789E+00  9.8962E-01  1.1360E+00  5.3463E+00  1.6958E+00  8.0758E-01  9.2031E-01
             9.6239E-01
 PARAMETER:  1.4285E-01 -4.5749E+00  7.8349E-01  6.1813E-01  8.9566E-02  2.2752E-01  1.7764E+00  6.2813E-01 -1.1371E-01  1.6950E-02
             6.1668E-02
 GRADIENT:   4.3684E+00  0.0000E+00 -1.6697E-02 -3.9575E+01 -1.7547E-01  1.2939E+00  2.5928E-03 -4.9594E-02 -1.3790E-01  6.5420E-02
            -7.1984E-03

0ITERATION NO.:  185    OBJECTIVE VALUE:  -1695.89994377618        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     6569             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0438E+00  1.0000E-02  1.9811E+00  1.6796E+00  9.8975E-01  1.1360E+00  5.8885E+00  1.6970E+00  8.0786E-01  9.1974E-01
             9.6241E-01
 PARAMETER:  1.4283E-01 -4.5731E+00  7.8364E-01  6.1853E-01  8.9702E-02  2.2750E-01  1.8730E+00  6.2887E-01 -1.1337E-01  1.6334E-02
             6.1683E-02
 GRADIENT:   7.5495E+02  0.0000E+00  1.2667E+01  1.5981E+03  6.1819E+00  1.7610E+02  2.6605E-01  3.6581E+00  2.1943E+01  6.1361E-01
             7.9151E-01

0ITERATION NO.:  190    OBJECTIVE VALUE:  -1695.90002050023        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     6752
 NPARAMETR:  1.0438E+00  1.0000E-02  1.9822E+00  1.6796E+00  9.9000E-01  1.1360E+00  5.8890E+00  1.6982E+00  8.0780E-01  9.1964E-01
             9.6241E-01
 PARAMETER:  1.4284E-01 -4.5731E+00  7.8421E-01  6.1856E-01  8.9953E-02  2.2751E-01  1.8731E+00  6.2956E-01 -1.1344E-01  1.6223E-02
             6.1685E-02
 GRADIENT:   4.2284E+00  0.0000E+00 -1.4287E-01 -3.8293E+01 -7.4149E-02  1.2483E+00  4.0079E-03  8.0303E-03 -1.7813E-02 -2.1771E-02
            -1.1633E-02

0ITERATION NO.:  195    OBJECTIVE VALUE:  -1695.90009161984        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     6929
 NPARAMETR:  1.0438E+00  1.0000E-02  1.9852E+00  1.6797E+00  9.9021E-01  1.1360E+00  6.2596E+00  1.6998E+00  8.0802E-01  9.2049E-01
             9.6249E-01
 PARAMETER:  1.4284E-01 -4.5731E+00  7.8470E-01  6.1859E-01  9.0249E-02  2.2751E-01  1.9150E+00  6.3007E-01 -1.1342E-01  1.6345E-02
             6.1678E-02
 GRADIENT:  -3.4245E-03  0.0000E+00 -1.0991E-01 -1.4166E-02  3.5411E-02  3.4089E-04 -6.3853E-05 -1.8124E-02 -2.9972E-02 -3.0473E-02
            -1.2748E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     6929
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.9765E-04 -1.4312E-03 -3.4615E-02 -6.2591E-03 -4.6039E-02
 SE:             2.9894E-02  9.4056E-04  1.9368E-02  2.9439E-02  1.9232E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9472E-01  1.2810E-01  7.3900E-02  8.3163E-01  1.6669E-02

 ETASHRINKSD(%)  1.0000E-10  9.6849E+01  3.5115E+01  1.3767E+00  3.5572E+01
 ETASHRINKVR(%)  1.0000E-10  9.9901E+01  5.7900E+01  2.7345E+00  5.8490E+01
 EBVSHRINKSD(%)  2.9575E-01  9.6979E+01  3.7659E+01  1.9106E+00  3.2809E+01
 EBVSHRINKVR(%)  5.9063E-01  9.9909E+01  6.1136E+01  3.7847E+00  5.4853E+01
 RELATIVEINF(%)  9.6807E+01  4.0306E-03  1.0708E+01  5.0588E+00  7.7530E+00
 EPSSHRINKSD(%)  4.5386E+01
 EPSSHRINKVR(%)  7.0173E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1695.9000916198361     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -960.74926505609790     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    33.68
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1695.900       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.00E-02  1.98E+00  1.68E+00  9.90E-01  1.14E+00  6.14E+00  1.70E+00  8.08E-01  9.20E-01  9.62E-01
 


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
 #CPUT: Total CPU Time in Seconds,      220.592
Stop Time:
Sun Oct 24 04:22:04 CDT 2021
