Sat Sep 18 10:12:07 CDT 2021
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
$DATA ../../../../data/spa/A2/dat100.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       18 SEP 2021
Days until program expires : 211
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
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -954.280993232289        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5612E+02 -1.7269E+01  3.5496E+01 -4.6893E+01  6.7073E+01 -4.2293E+00 -2.1446E+01 -2.5834E+01 -3.1857E+01 -6.3168E+01
            -1.2126E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1338.92733251880        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.8910E-01  1.0039E+00  9.5634E-01  1.0418E+00  9.3835E-01  9.9291E-01  1.0454E+00  1.0349E+00  1.0326E+00  1.0778E+00
             1.8404E+00
 PARAMETER:  8.9037E-02  1.0391E-01  5.5358E-02  1.4094E-01  3.6363E-02  9.2885E-02  1.4442E-01  1.3431E-01  1.3212E-01  1.7491E-01
             7.0998E-01
 GRADIENT:   8.5564E+01  6.0881E+00  9.7238E+00  1.8123E-01 -1.8156E+00 -2.2554E+00 -3.4054E+00 -3.7474E+00 -2.2737E+00 -5.9755E+00
            -2.0819E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1364.15807816195        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.4898E-01  8.4023E-01  1.2450E+00  1.1691E+00  9.8837E-01  9.5816E-01  1.3219E+00  1.2630E+00  9.2679E-01  1.0623E+00
             2.2034E+00
 PARAMETER:  4.7636E-02 -7.4078E-02  3.1916E-01  2.5621E-01  8.8298E-02  5.7255E-02  3.7907E-01  3.3348E-01  2.3968E-02  1.6047E-01
             8.9001E-01
 GRADIENT:  -1.2838E+01  1.2539E+01  9.6747E+00  1.2428E+01 -7.2870E+00 -1.0092E+01  4.9629E+00 -2.5074E+00  2.5496E+00 -5.7156E+00
            -7.5264E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1369.85057450051        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.6092E-01  8.7987E-01  1.0502E+00  1.1316E+00  9.3759E-01  9.8107E-01  1.1086E+00  9.3791E-01  9.1811E-01  1.0211E+00
             2.5415E+00
 PARAMETER:  6.0134E-02 -2.7976E-02  1.4901E-01  2.2365E-01  3.5557E-02  8.0893E-02  2.0311E-01  3.5895E-02  1.4558E-02  1.2089E-01
             1.0327E+00
 GRADIENT:   1.8756E+00 -1.7251E+00 -3.9707E+00 -5.5491E-01  2.1871E+00  3.8708E-01 -3.3702E-01  1.7772E+00 -9.2688E-01  5.3244E-01
             4.2186E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1371.45243676505        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  9.5586E-01  6.8397E-01  1.1806E+00  1.2591E+00  9.2564E-01  9.7214E-01  1.0665E+00  1.4066E-01  9.0215E-01  1.0737E+00
             2.5704E+00
 PARAMETER:  5.4855E-02 -2.7983E-01  2.6605E-01  3.3043E-01  2.2733E-02  7.1746E-02  1.6440E-01 -1.8614E+00 -2.9750E-03  1.7114E-01
             1.0441E+00
 GRADIENT:  -3.8846E+00  1.0488E+00 -1.7488E+00  4.5631E+00  2.0373E+00 -9.5664E-01 -9.3012E-01  6.7120E-02  1.5664E-01 -8.6495E-01
             7.4454E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1371.91964103158        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  9.5396E-01  4.7501E-01  1.1717E+00  1.3854E+00  8.4819E-01  9.7027E-01  1.4302E+00  1.3302E-02  8.3217E-01  1.0570E+00
             2.5591E+00
 PARAMETER:  5.2870E-02 -6.4441E-01  2.5843E-01  4.2598E-01 -6.4653E-02  6.9818E-02  4.5782E-01 -4.2198E+00 -8.3716E-02  1.5546E-01
             1.0396E+00
 GRADIENT:  -3.1441E+00  2.4349E+00  7.7767E-01  1.0082E+01 -2.7595E+00 -9.7151E-01 -2.5852E-01  6.6905E-04 -1.9043E-01  2.5623E-01
             3.5125E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1372.08602311153        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      552
 NPARAMETR:  9.5639E-01  3.6550E-01  1.2356E+00  1.4582E+00  8.4528E-01  9.7343E-01  1.8628E+00  1.0000E-02  7.9516E-01  1.0600E+00
             2.5668E+00
 PARAMETER:  5.5410E-02 -9.0648E-01  3.1154E-01  4.7717E-01 -6.8087E-02  7.3069E-02  7.2209E-01 -6.0289E+00 -1.2921E-01  1.5824E-01
             1.0427E+00
 GRADIENT:   3.9447E-01  2.1872E+00  1.0021E+00  7.1272E+00 -1.4740E+00  2.5902E-01  7.7504E-01  0.0000E+00  6.7733E-02 -4.6330E-01
            -3.1013E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1372.29211599988        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      728
 NPARAMETR:  9.5373E-01  2.1676E-01  1.2439E+00  1.5416E+00  8.0922E-01  9.7123E-01  2.1268E+00  1.0000E-02  7.7354E-01  1.0625E+00
             2.5639E+00
 PARAMETER:  5.2623E-02 -1.4289E+00  3.1822E-01  5.3284E-01 -1.1169E-01  7.0804E-02  8.5460E-01 -1.0980E+01 -1.5678E-01  1.6064E-01
             1.0415E+00
 GRADIENT:  -1.3507E-01  2.2927E-01 -1.2453E-01  9.8965E-01 -9.4615E-02 -1.5435E-02 -1.3928E-01  0.0000E+00 -4.9229E-01  1.9449E-02
            -3.6636E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1372.34804565499        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      903
 NPARAMETR:  9.5237E-01  1.2407E-01  1.2629E+00  1.5980E+00  7.9371E-01  9.7009E-01  2.6262E+00  1.0000E-02  7.5604E-01  1.0620E+00
             2.5644E+00
 PARAMETER:  5.1201E-02 -1.9869E+00  3.3339E-01  5.6873E-01 -1.3104E-01  6.9630E-02  1.0655E+00 -1.6425E+01 -1.7966E-01  1.6017E-01
             1.0417E+00
 GRADIENT:   3.1148E-01  1.5700E-01 -5.8055E-01  3.1370E+00  4.0413E-01 -6.5440E-02 -1.4457E-01  0.0000E+00 -9.7779E-01 -1.8469E-01
            -7.2340E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1372.38608638965        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1080
 NPARAMETR:  9.5107E-01  6.4258E-02  1.2850E+00  1.6322E+00  7.8740E-01  9.6896E-01  4.1083E+00  1.0000E-02  7.4573E-01  1.0645E+00
             2.5678E+00
 PARAMETER:  4.9834E-02 -2.6448E+00  3.5074E-01  5.8990E-01 -1.3901E-01  6.8469E-02  1.5130E+00 -2.2692E+01 -1.9339E-01  1.6247E-01
             1.0430E+00
 GRADIENT:  -3.5744E-02 -1.1119E-02 -9.9129E-02 -1.2416E+00  3.9196E-01 -3.2719E-02  1.5644E-03  0.0000E+00  2.0964E-03 -6.5263E-02
             2.1008E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1372.39091219836        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1255
 NPARAMETR:  9.5066E-01  4.2159E-02  1.2829E+00  1.6456E+00  7.8067E-01  9.6876E-01  5.4901E+00  1.0000E-02  7.4221E-01  1.0645E+00
             2.5663E+00
 PARAMETER:  4.9396E-02 -3.0663E+00  3.4916E-01  5.9808E-01 -1.4760E-01  6.8263E-02  1.8030E+00 -2.6775E+01 -1.9812E-01  1.6251E-01
             1.0425E+00
 GRADIENT:  -1.2599E-01  6.4514E-02  2.3998E-01  6.8176E-01 -5.2402E-01 -1.9120E-02  6.8278E-02  0.0000E+00  8.7120E-02  9.2731E-02
             1.0988E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1372.39897331936        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1431
 NPARAMETR:  9.5022E-01  1.7792E-02  1.2960E+00  1.6604E+00  7.8021E-01  9.6839E-01  7.4987E+00  1.0000E-02  7.3813E-01  1.0670E+00
             2.5670E+00
 PARAMETER:  4.8943E-02 -3.9290E+00  3.5929E-01  6.0705E-01 -1.4819E-01  6.7876E-02  2.1147E+00 -3.5799E+01 -2.0363E-01  1.6486E-01
             1.0427E+00
 GRADIENT:  -7.8272E-02  9.3211E-04  4.6829E-02 -4.5647E-02 -8.4749E-02 -3.9125E-02 -7.0741E-03  0.0000E+00 -6.3912E-02  4.9834E-02
            -1.1215E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1372.39976356116        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1607
 NPARAMETR:  9.5013E-01  1.0000E-02  1.2975E+00  1.6654E+00  7.7885E-01  9.6845E-01  9.0982E+00  1.0000E-02  7.3697E-01  1.0663E+00
             2.5671E+00
 PARAMETER:  4.8844E-02 -4.5178E+00  3.6041E-01  6.1008E-01 -1.4993E-01  6.7940E-02  2.3081E+00 -4.2037E+01 -2.0520E-01  1.6423E-01
             1.0428E+00
 GRADIENT:   3.3944E-03  0.0000E+00  4.4777E-02  8.0639E-01 -1.7182E-01  8.7237E-03 -9.0252E-03  0.0000E+00 -7.3862E-02 -1.4553E-02
            -6.4771E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1372.40061117531        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1783
 NPARAMETR:  9.5011E-01  1.0000E-02  1.2972E+00  1.6650E+00  7.7890E-01  9.6838E-01  1.0419E+01  1.0000E-02  7.3697E-01  1.0665E+00
             2.5671E+00
 PARAMETER:  4.8821E-02 -4.8463E+00  3.6024E-01  6.0983E-01 -1.4988E-01  6.7868E-02  2.4437E+00 -4.5486E+01 -2.0521E-01  1.6435E-01
             1.0428E+00
 GRADIENT:  -5.5396E-03  0.0000E+00 -7.0632E-05 -4.3094E-03  1.9092E-03 -1.3602E-03  1.0716E-04  0.0000E+00  1.4535E-03  1.0331E-03
             2.4371E-04

0ITERATION NO.:   69    OBJECTIVE VALUE:  -1372.40061252312        NO. OF FUNC. EVALS.: 129
 CUMULATIVE NO. OF FUNC. EVALS.:     1912
 NPARAMETR:  9.5011E-01  1.0000E-02  1.2972E+00  1.6650E+00  7.7889E-01  9.6839E-01  1.0412E+01  1.0000E-02  7.3697E-01  1.0665E+00
             2.5671E+00
 PARAMETER:  4.8824E-02 -4.8463E+00  3.6024E-01  6.0983E-01 -1.4989E-01  6.7876E-02  2.4429E+00 -4.5486E+01 -2.0521E-01  1.6439E-01
             1.0428E+00
 GRADIENT:  -3.4377E-04  0.0000E+00  3.8235E-03 -7.5936E-03 -6.6922E-03  1.5144E-03  1.9631E-05  0.0000E+00  1.7351E-03  6.1064E-03
             5.8756E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1912
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.8808E-05 -3.7292E-04  6.6798E-05 -1.5184E-02 -3.0508E-02
 SE:             2.9050E-02  1.5520E-03  1.0680E-04  2.6396E-02  1.9991E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9838E-01  8.1011E-01  5.3168E-01  5.6513E-01  1.2698E-01

 ETASHRINKSD(%)  2.6802E+00  9.4801E+01  9.9642E+01  1.1569E+01  3.3029E+01
 ETASHRINKVR(%)  5.2885E+00  9.9730E+01  9.9999E+01  2.1799E+01  5.5148E+01
 EBVSHRINKSD(%)  2.4745E+00  9.4903E+01  9.9616E+01  1.0955E+01  3.3409E+01
 EBVSHRINKVR(%)  4.8878E+00  9.9740E+01  9.9999E+01  2.0711E+01  5.5656E+01
 RELATIVEINF(%)  8.2833E+01  2.6686E-03  9.3091E-05  1.0284E+00  1.8896E+00
 EPSSHRINKSD(%)  3.1039E+01
 EPSSHRINKVR(%)  5.2444E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1372.4006125231160     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -637.24978595937785     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.87
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.27
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1372.401       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.50E-01  1.00E-02  1.30E+00  1.67E+00  7.79E-01  9.68E-01  1.04E+01  1.00E-02  7.37E-01  1.07E+00  2.57E+00
 


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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.30E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.45E+00  0.00E+00  9.85E+01
 
 TH 4
+       -4.92E+01  0.00E+00 -1.79E+01  5.87E+02
 
 TH 5
+        1.77E+01  0.00E+00 -2.43E+02 -9.53E+01  6.95E+02
 
 TH 6
+        1.30E+01  0.00E+00  1.06E+01 -1.57E+01 -5.56E+00  1.99E+02
 
 TH 7
+        5.84E-03  0.00E+00  6.80E-03 -4.78E-03  9.18E-03 -1.48E-02  1.08E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.57E+01  0.00E+00  1.46E+01 -1.06E+01 -2.77E+00  6.06E+00  9.64E-02  0.00E+00  2.39E+02
 
 TH10
+       -1.32E+01  0.00E+00  3.92E+00 -2.29E+00 -3.05E+01 -1.43E+00  2.13E-02  0.00E+00 -1.43E+00  5.30E+01
 
 TH11
+       -1.47E+01  0.00E+00 -4.03E+00 -1.21E+01 -5.72E+00  2.05E+00  3.24E-03  0.00E+00  1.37E+01  1.45E+01  4.62E+01
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       29.205
Stop Time:
Sat Sep 18 10:12:38 CDT 2021
