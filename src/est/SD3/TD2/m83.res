Sun Oct 24 00:50:07 CDT 2021
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
$DATA ../../../../data/SD3/TD2/dat83.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m83.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2132.22920484923        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3239E+02  2.7591E+01 -2.3472E-01  4.1578E+01 -3.1443E+01  2.7287E+01  1.3110E+01  9.9666E+00  1.5298E+01  2.6095E+01
            -6.7934E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2135.42570242642        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      196
 NPARAMETR:  1.0074E+00  9.9321E-01  1.1201E+00  1.0260E+00  1.0779E+00  1.0730E+00  9.1959E-01  9.2890E-01  9.5527E-01  8.3580E-01
             1.0210E+00
 PARAMETER:  1.0740E-01  9.3187E-02  2.1340E-01  1.2562E-01  1.7503E-01  1.7046E-01  1.6168E-02  2.6243E-02  5.4236E-02 -7.9372E-02
             1.2078E-01
 GRADIENT:   5.2507E+00 -2.7437E+00  8.1675E+00 -1.4239E+01  2.3870E+01  1.1718E+01  8.8118E-02 -5.1563E+00 -7.1566E+00 -1.2763E+01
            -1.4358E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2136.11743497440        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  9.9810E-01  8.3781E-01  1.1593E+00  1.1388E+00  1.0082E+00  1.0799E+00  7.9870E-01  8.7524E-01  9.3185E-01  8.4037E-01
             1.0234E+00
 PARAMETER:  9.8095E-02 -7.6961E-02  2.4785E-01  2.2999E-01  1.0816E-01  1.7684E-01 -1.2477E-01 -3.3258E-02  2.9412E-02 -7.3910E-02
             1.2316E-01
 GRADIENT:  -1.0029E+01  1.9476E+01  1.7801E+01  1.7858E+01 -7.2652E+00  1.4513E+01 -2.7307E+00 -7.3793E+00 -6.6189E+00 -9.9014E+00
             2.0217E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2137.52006013210        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  1.0021E+00  7.2289E-01  1.1814E+00  1.2041E+00  9.7382E-01  1.0391E+00  9.2611E-01  9.3172E-01  8.9666E-01  8.8193E-01
             1.0162E+00
 PARAMETER:  1.0213E-01 -2.2450E-01  2.6671E-01  2.8575E-01  7.3471E-02  1.3836E-01  2.3237E-02  2.9279E-02 -9.0790E-03 -2.5641E-02
             1.1610E-01
 GRADIENT:   5.2777E-01  7.8427E+00  6.6957E+00  6.5437E+00 -1.2167E+01  2.1701E-01  2.1443E-01 -7.0636E-01 -3.0876E-01  1.6614E+00
             5.9261E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2137.73644037159        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      723
 NPARAMETR:  9.9924E-01  5.6418E-01  1.2802E+00  1.3056E+00  9.6484E-01  1.0379E+00  9.0859E-01  1.0037E+00  8.5166E-01  8.8286E-01
             1.0152E+00
 PARAMETER:  9.9236E-02 -4.7239E-01  3.4704E-01  3.6669E-01  6.4206E-02  1.3717E-01  4.1402E-03  1.0365E-01 -6.0567E-02 -2.4593E-02
             1.1510E-01
 GRADIENT:  -4.4863E-01  4.2828E+00  2.0757E+00  8.9033E+00 -2.9021E+00  4.4908E-01 -8.3262E-02 -3.3442E-01 -3.0299E-01 -1.1417E-01
            -3.1909E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2137.77404393675        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      905
 NPARAMETR:  9.9905E-01  5.1173E-01  1.2997E+00  1.3316E+00  9.5863E-01  1.0364E+00  9.2644E-01  1.0205E+00  8.3808E-01  8.8334E-01
             1.0153E+00
 PARAMETER:  9.9049E-02 -5.6995E-01  3.6212E-01  3.8635E-01  5.7745E-02  1.3579E-01  2.3595E-02  1.2029E-01 -7.6638E-02 -2.4047E-02
             1.1522E-01
 GRADIENT:   9.9755E-01 -7.5861E-01 -7.8428E-01 -5.5811E+00  2.1816E+00  1.5520E-01  6.7907E-02  1.1881E-01  7.3052E-01 -8.7529E-02
             6.4626E-04

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2137.77684328995        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1089             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9936E-01  5.1164E-01  1.2995E+00  1.3321E+00  9.5740E-01  1.0365E+00  9.2538E-01  1.0193E+00  8.3654E-01  8.8325E-01
             1.0154E+00
 PARAMETER:  9.9355E-02 -5.7013E-01  3.6196E-01  3.8673E-01  5.6468E-02  1.3587E-01  2.2450E-02  1.1915E-01 -7.8480E-02 -2.4151E-02
             1.1527E-01
 GRADIENT:   4.2823E+02  6.2215E+01  7.6450E+00  5.1966E+02  6.7171E+00  5.6111E+01  8.0472E-01  3.3308E-01  1.0900E+01  5.9775E-01
             1.1576E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2137.77714657170        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1265
 NPARAMETR:  9.9940E-01  5.1182E-01  1.2987E+00  1.3323E+00  9.5688E-01  1.0365E+00  9.2278E-01  1.0185E+00  8.3616E-01  8.8274E-01
             1.0154E+00
 PARAMETER:  9.9397E-02 -5.6979E-01  3.6136E-01  3.8687E-01  5.5923E-02  1.3589E-01  1.9635E-02  1.1835E-01 -7.8937E-02 -2.4721E-02
             1.1527E-01
 GRADIENT:   1.6782E+00  1.3762E-01  3.0414E-01 -3.8313E+00  3.9643E-02  1.9277E-01 -1.7899E-02  9.3493E-03 -9.3628E-02  1.7913E-02
            -1.8468E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2137.77741729735        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1444
 NPARAMETR:  9.9935E-01  5.1183E-01  1.2979E+00  1.3321E+00  9.5652E-01  1.0365E+00  9.2466E-01  1.0177E+00  8.3608E-01  8.8243E-01
             1.0154E+00
 PARAMETER:  9.9353E-02 -5.6976E-01  3.6078E-01  3.8678E-01  5.5551E-02  1.3586E-01  2.1668E-02  1.1755E-01 -7.9034E-02 -2.5075E-02
             1.1528E-01
 GRADIENT:   1.5883E+00  1.0577E-01  3.3795E-01 -4.0011E+00 -2.1924E-02  1.8303E-01 -1.2698E-02 -3.7663E-03 -8.8372E-02  1.2938E-02
            -2.0376E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2137.77761224932        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1629
 NPARAMETR:  9.9936E-01  5.1184E-01  1.2972E+00  1.3320E+00  9.5621E-01  1.0365E+00  9.2363E-01  1.0169E+00  8.3611E-01  8.8219E-01
             1.0154E+00
 PARAMETER:  9.9356E-02 -5.6975E-01  3.6022E-01  3.8671E-01  5.5225E-02  1.3587E-01  2.0554E-02  1.1680E-01 -7.8993E-02 -2.5353E-02
             1.1528E-01
 GRADIENT:   1.5903E+00  8.4603E-02  3.2874E-01 -4.0865E+00 -3.5690E-02  1.8363E-01 -1.9276E-02 -1.3067E-02 -1.1014E-01  5.8972E-04
            -1.8683E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2137.77782391060        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1812
 NPARAMETR:  9.9936E-01  5.1190E-01  1.2965E+00  1.3320E+00  9.5594E-01  1.0365E+00  9.2628E-01  1.0164E+00  8.3607E-01  8.8194E-01
             1.0154E+00
 PARAMETER:  9.9358E-02 -5.6962E-01  3.5969E-01  3.8665E-01  5.4940E-02  1.3587E-01  2.3421E-02  1.1622E-01 -7.9047E-02 -2.5634E-02
             1.1528E-01
 GRADIENT:   1.5898E+00  7.3925E-02  3.0730E-01 -4.0978E+00 -3.9162E-02  1.8342E-01 -9.4448E-03 -1.1349E-02 -7.5118E-02  4.8799E-03
            -2.2496E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2137.77806611464        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1996             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9936E-01  5.1203E-01  1.2955E+00  1.3318E+00  9.5574E-01  1.0365E+00  9.2866E-01  1.0158E+00  8.3632E-01  8.8176E-01
             1.0154E+00
 PARAMETER:  9.9364E-02 -5.6936E-01  3.5892E-01  3.8650E-01  5.4735E-02  1.3587E-01  2.5986E-02  1.1571E-01 -7.8744E-02 -2.5836E-02
             1.1530E-01
 GRADIENT:   4.2812E+02  6.2275E+01  7.5580E+00  5.1939E+02  6.5559E+00  5.6094E+01  8.0416E-01  3.1348E-01  1.0815E+01  5.7077E-01
             1.1577E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2137.77818023229        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2179
 NPARAMETR:  9.9937E-01  5.1217E-01  1.2949E+00  1.3316E+00  9.5557E-01  1.0365E+00  9.2982E-01  1.0154E+00  8.3643E-01  8.8159E-01
             1.0154E+00
 PARAMETER:  9.9369E-02 -5.6909E-01  3.5842E-01  3.8637E-01  5.4552E-02  1.3587E-01  2.7238E-02  1.1528E-01 -7.8609E-02 -2.6026E-02
             1.1530E-01
 GRADIENT:   1.6040E+00 -9.5122E-02  1.7857E-02 -4.3799E+00  3.0824E-01  1.8273E-01  1.4795E-02  3.6949E-02  1.0547E-01  1.6862E-02
             1.0871E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2137.77836450196        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2363             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9937E-01  5.1258E-01  1.2943E+00  1.3315E+00  9.5529E-01  1.0365E+00  9.2688E-01  1.0147E+00  8.3633E-01  8.8138E-01
             1.0154E+00
 PARAMETER:  9.9371E-02 -5.6830E-01  3.5800E-01  3.8629E-01  5.4259E-02  1.3587E-01  2.4069E-02  1.1455E-01 -7.8726E-02 -2.6267E-02
             1.1529E-01
 GRADIENT:   4.2810E+02  6.2390E+01  7.6235E+00  5.1902E+02  6.3489E+00  5.6095E+01  7.8586E-01  2.9328E-01  1.0708E+01  5.6093E-01
             1.1431E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2137.77841063061        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2542
 NPARAMETR:  9.9937E-01  5.1285E-01  1.2939E+00  1.3314E+00  9.5511E-01  1.0365E+00  9.2563E-01  1.0141E+00  8.3631E-01  8.8123E-01
             1.0154E+00
 PARAMETER:  9.9373E-02 -5.6777E-01  3.5765E-01  3.8620E-01  5.4067E-02  1.3587E-01  2.2714E-02  1.1404E-01 -7.8762E-02 -2.6439E-02
             1.1529E-01
 GRADIENT:   1.5809E+00  1.4202E-01  2.2418E-01 -3.8051E+00 -1.1411E-01  1.8132E-01 -2.1063E-02 -5.9972E-03 -1.2597E-01 -8.4440E-04
            -1.9640E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2137.77858862398        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2726             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9938E-01  5.1289E-01  1.2932E+00  1.3311E+00  9.5502E-01  1.0365E+00  9.3057E-01  1.0138E+00  8.3661E-01  8.8112E-01
             1.0154E+00
 PARAMETER:  9.9382E-02 -5.6770E-01  3.5708E-01  3.8599E-01  5.3979E-02  1.3588E-01  2.8044E-02  1.1371E-01 -7.8401E-02 -2.6557E-02
             1.1531E-01
 GRADIENT:   4.2809E+02  6.2245E+01  7.4686E+00  5.1815E+02  6.5464E+00  5.6093E+01  8.1586E-01  3.1224E-01  1.0850E+01  5.7145E-01
             1.1628E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2137.77868253770        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2909
 NPARAMETR:  9.9938E-01  5.1308E-01  1.2929E+00  1.3310E+00  9.5491E-01  1.0365E+00  9.2938E-01  1.0135E+00  8.3661E-01  8.8104E-01
             1.0154E+00
 PARAMETER:  9.9384E-02 -5.6732E-01  3.5686E-01  3.8593E-01  5.3864E-02  1.3588E-01  2.6766E-02  1.1338E-01 -7.8401E-02 -2.6657E-02
             1.1530E-01
 GRADIENT:   1.5975E+00 -2.7462E-02  5.2820E-02 -4.2042E+00  1.3151E-01  1.8244E-01  4.9152E-03  2.0450E-02  4.5051E-02  1.0799E-02
             3.4885E-03

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2137.77877907733        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3093             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9939E-01  5.1350E-01  1.2924E+00  1.3308E+00  9.5470E-01  1.0365E+00  9.2764E-01  1.0128E+00  8.3657E-01  8.8087E-01
             1.0154E+00
 PARAMETER:  9.9388E-02 -5.6650E-01  3.5648E-01  3.8580E-01  5.3641E-02  1.3588E-01  2.4891E-02  1.1274E-01 -7.8441E-02 -2.6845E-02
             1.1530E-01
 GRADIENT:   4.2807E+02  6.2424E+01  7.6160E+00  5.1792E+02  6.2392E+00  5.6096E+01  7.8874E-01  2.7756E-01  1.0694E+01  5.6020E-01
             1.1410E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2137.77886093047        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     3274
 NPARAMETR:  9.9939E-01  5.1359E-01  1.2920E+00  1.3307E+00  9.5466E-01  1.0365E+00  9.2867E-01  1.0126E+00  8.3667E-01  8.8081E-01
             1.0154E+00
 PARAMETER:  9.9392E-02 -5.6632E-01  3.5622E-01  3.8570E-01  5.3596E-02  1.3588E-01  2.5993E-02  1.1255E-01 -7.8325E-02 -2.6916E-02
             1.1530E-01
 GRADIENT:   1.5903E+00  3.2453E-02  9.9025E-02 -4.0563E+00  8.9508E-03  1.8200E-01 -3.4000E-03  7.7131E-03 -1.0231E-02  6.0120E-03
            -4.3646E-03

0ITERATION NO.:   95    OBJECTIVE VALUE:  -2137.77894157888        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3458             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9940E-01  5.1378E-01  1.2914E+00  1.3304E+00  9.5459E-01  1.0365E+00  9.3052E-01  1.0123E+00  8.3684E-01  8.8069E-01
             1.0154E+00
 PARAMETER:  9.9400E-02 -5.6595E-01  3.5574E-01  3.8551E-01  5.3523E-02  1.3589E-01  2.7987E-02  1.1220E-01 -7.8121E-02 -2.7052E-02
             1.1531E-01
 GRADIENT:   4.2807E+02  6.2252E+01  7.4047E+00  5.1699E+02  6.5462E+00  5.6096E+01  8.1440E-01  3.0527E-01  1.0819E+01  5.6212E-01
             1.1619E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -2137.77900881249        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     3641
 NPARAMETR:  9.9940E-01  5.1399E-01  1.2912E+00  1.3304E+00  9.5450E-01  1.0365E+00  9.3016E-01  1.0120E+00  8.3684E-01  8.8063E-01
             1.0154E+00
 PARAMETER:  9.9401E-02 -5.6554E-01  3.5560E-01  3.8545E-01  5.3436E-02  1.3589E-01  2.7604E-02  1.1194E-01 -7.8123E-02 -2.7116E-02
             1.1531E-01
 GRADIENT:   1.5980E+00 -2.3637E-02  2.7891E-02 -4.1807E+00  9.4437E-02  1.8269E-01  5.2416E-03  1.5997E-02  4.0023E-02  7.5081E-03
             2.2272E-03

0ITERATION NO.:  105    OBJECTIVE VALUE:  -2137.77908223328        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3825             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9940E-01  5.1439E-01  1.2909E+00  1.3302E+00  9.5436E-01  1.0365E+00  9.2883E-01  1.0115E+00  8.3686E-01  8.8053E-01
             1.0154E+00
 PARAMETER:  9.9405E-02 -5.6477E-01  3.5532E-01  3.8533E-01  5.3284E-02  1.3589E-01  2.6175E-02  1.1145E-01 -7.8098E-02 -2.7228E-02
             1.1530E-01
 GRADIENT:   4.2806E+02  6.2414E+01  7.5719E+00  5.1677E+02  6.2298E+00  5.6099E+01  7.9824E-01  2.7299E-01  1.0718E+01  5.6160E-01
             1.1440E+00

0ITERATION NO.:  110    OBJECTIVE VALUE:  -2137.77913788286        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     4006
 NPARAMETR:  9.9941E-01  5.1449E-01  1.2906E+00  1.3301E+00  9.5434E-01  1.0365E+00  9.2946E-01  1.0114E+00  8.3694E-01  8.8049E-01
             1.0154E+00
 PARAMETER:  9.9408E-02 -5.6459E-01  3.5512E-01  3.8523E-01  5.3268E-02  1.3590E-01  2.6846E-02  1.1132E-01 -7.8006E-02 -2.7281E-02
             1.1530E-01
 GRADIENT:   1.5931E+00  2.3363E-02  7.1552E-02 -4.0694E+00 -9.0344E-03  1.8271E-01 -1.0334E-03  6.0054E-03  1.7877E-03  5.3183E-03
            -2.9532E-03

0ITERATION NO.:  111    OBJECTIVE VALUE:  -2137.77913788286        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:     4030
 NPARAMETR:  9.9941E-01  5.1449E-01  1.2906E+00  1.3301E+00  9.5434E-01  1.0365E+00  9.2946E-01  1.0114E+00  8.3694E-01  8.8049E-01
             1.0154E+00
 PARAMETER:  9.9408E-02 -5.6459E-01  3.5512E-01  3.8523E-01  5.3268E-02  1.3590E-01  2.6846E-02  1.1132E-01 -7.8006E-02 -2.7281E-02
             1.1530E-01
 GRADIENT:  -5.6130E-03 -2.0407E-02  9.4385E-02  1.9179E-01 -1.4954E-04 -6.4786E-04 -1.9595E-03  5.6734E-03 -2.4232E-02  4.6513E-03
            -2.4230E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     4030
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.9832E-04 -1.1315E-02 -2.8128E-02 -4.9218E-04 -3.2275E-02
 SE:             2.9883E-02  8.6089E-03  1.6788E-02  2.8461E-02  2.1246E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8936E-01  1.8872E-01  9.3839E-02  9.8620E-01  1.2873E-01

 ETASHRINKSD(%)  1.0000E-10  7.1159E+01  4.3757E+01  4.6520E+00  2.8825E+01
 ETASHRINKVR(%)  1.0000E-10  9.1682E+01  6.8367E+01  9.0875E+00  4.9341E+01
 EBVSHRINKSD(%)  3.1484E-01  7.1771E+01  4.6359E+01  4.8158E+00  2.7122E+01
 EBVSHRINKVR(%)  6.2869E-01  9.2031E+01  7.1226E+01  9.3997E+00  4.6887E+01
 RELATIVEINF(%)  9.7460E+01  3.2627E-01  6.1885E+00  4.6540E+00  7.9167E+00
 EPSSHRINKSD(%)  3.3581E+01
 EPSSHRINKVR(%)  5.5885E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2137.7791378828592     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1218.8406046781865     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.46
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2137.779       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  5.14E-01  1.29E+00  1.33E+00  9.54E-01  1.04E+00  9.29E-01  1.01E+00  8.37E-01  8.80E-01  1.02E+00
 


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
 #CPUT: Total CPU Time in Seconds,      142.886
Stop Time:
Sun Oct 24 00:50:31 CDT 2021
