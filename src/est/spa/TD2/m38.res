Wed Sep 29 18:58:41 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat38.csv ignore=@
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
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
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
 RAW OUTPUT FILE (FILE): m38.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1614.63551224275        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6495E+02  1.9606E+01 -1.5967E+01  5.1290E+01 -1.3021E+01  3.1054E+01 -2.1008E+01  1.2558E+01 -3.7810E+01  6.1453E+00
            -3.6117E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1625.62703710971        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  1.0025E+00  1.0536E+00  1.1385E+00  9.8976E-01  1.0740E+00  1.0831E+00  1.1239E+00  9.3050E-01  1.2003E+00  9.8207E-01
             1.1034E+00
 PARAMETER:  1.0247E-01  1.5218E-01  2.2969E-01  8.9704E-02  1.7139E-01  1.7987E-01  2.1679E-01  2.7962E-02  2.8260E-01  8.1908E-02
             1.9837E-01
 GRADIENT:   1.8492E+01  3.1883E+00  5.6853E+00  2.9927E-02 -2.9087E+01  6.2988E+00 -3.2670E+00  5.0344E+00  1.8887E+00 -3.3449E+00
             6.0219E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1626.57246237746        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.9431E-01  1.1453E+00  1.1748E+00  9.4591E-01  1.1953E+00  1.0307E+00  1.0947E+00  7.2318E-01  1.2865E+00  1.1635E+00
             1.0629E+00
 PARAMETER:  9.4298E-02  2.3563E-01  2.6112E-01  4.4393E-02  2.7842E-01  1.3025E-01  1.9050E-01 -2.2409E-01  3.5190E-01  2.5147E-01
             1.6104E-01
 GRADIENT:   3.4100E+00  1.0144E+00 -2.3934E+00  1.7954E+01  1.7028E+01 -1.3978E+01  2.1258E+00  5.3101E-02  5.8835E+00  3.0891E+00
            -8.7489E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1628.80034951568        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  9.9653E-01  1.3101E+00  9.1069E-01  8.2591E-01  1.1298E+00  1.0698E+00  1.0485E+00  3.5656E-01  1.3241E+00  1.0350E+00
             1.0840E+00
 PARAMETER:  9.6521E-02  3.7012E-01  6.4459E-03 -9.1271E-02  2.2204E-01  1.6747E-01  1.4739E-01 -9.3125E-01  3.8076E-01  1.3440E-01
             1.8064E-01
 GRADIENT:   3.6427E+00  7.8876E+00 -9.1015E-01  9.2560E+00 -6.5257E+00  6.9933E-01 -1.3445E-01  6.9423E-01 -1.9346E+00  6.6513E-01
             4.2482E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1629.77883771201        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  9.9461E-01  1.6755E+00  7.9462E-01  5.9501E-01  1.2988E+00  1.0706E+00  8.5888E-01  1.5706E-01  1.7346E+00  1.1322E+00
             1.0906E+00
 PARAMETER:  9.4594E-02  6.1613E-01 -1.2989E-01 -4.1917E-01  3.6146E-01  1.6821E-01 -5.2131E-02 -1.7512E+00  6.5076E-01  2.2413E-01
             1.8676E-01
 GRADIENT:  -1.9273E+00  2.0166E+01  2.0674E+00  1.2829E+01 -1.2915E+00  7.6394E-01 -1.2024E+00  8.5263E-02 -9.3178E-01 -9.3406E-01
            -1.0467E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1630.08200817123        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      886
 NPARAMETR:  9.9579E-01  1.7429E+00  7.5134E-01  5.3297E-01  1.3306E+00  1.0687E+00  8.3950E-01  5.4485E-02  1.8547E+00  1.1521E+00
             1.0917E+00
 PARAMETER:  9.5785E-02  6.5557E-01 -1.8590E-01 -5.2928E-01  3.8565E-01  1.6642E-01 -7.4954E-02 -2.8098E+00  7.1772E-01  2.4161E-01
             1.8776E-01
 GRADIENT:   2.9798E-01 -1.1583E+00  1.0228E+00  1.9338E+00 -6.1963E-01  1.0029E-01  1.9461E-01  1.2106E-02 -3.6590E-01  4.6483E-02
            -4.6506E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1630.09773724765        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1068
 NPARAMETR:  9.9629E-01  1.7487E+00  7.4312E-01  5.2642E-01  1.3310E+00  1.0696E+00  8.3681E-01  1.4709E-02  1.8679E+00  1.1498E+00
             1.0925E+00
 PARAMETER:  9.6282E-02  6.5887E-01 -1.9689E-01 -5.4165E-01  3.8593E-01  1.6724E-01 -7.8153E-02 -4.1193E+00  7.2479E-01  2.3957E-01
             1.8847E-01
 GRADIENT:   1.2332E+00 -4.8749E+00  5.3444E-01  5.6531E-01 -8.0130E-01  4.3686E-01  1.4758E-01  9.3026E-04 -1.5950E-01  3.4467E-02
            -9.1242E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1630.10286902048        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1251
 NPARAMETR:  9.9648E-01  1.7461E+00  7.4325E-01  5.2565E-01  1.3303E+00  1.0699E+00  8.3587E-01  1.0000E-02  1.8706E+00  1.1488E+00
             1.0922E+00
 PARAMETER:  9.6472E-02  6.5741E-01 -1.9672E-01 -5.4312E-01  3.8537E-01  1.6761E-01 -7.9284E-02 -5.0187E+00  7.2623E-01  2.3868E-01
             1.8819E-01
 GRADIENT:   1.6519E+00 -8.1904E+00  4.1783E-01 -6.6774E-01 -7.1507E-01  5.8901E-01  6.9318E-02  0.0000E+00  6.3022E-02  1.6829E-02
            -1.5967E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1630.10434001144        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1445             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9646E-01  1.7441E+00  7.4180E-01  5.2680E-01  1.3301E+00  1.0699E+00  8.3585E-01  1.0000E-02  1.8699E+00  1.1481E+00
             1.0925E+00
 PARAMETER:  9.6457E-02  6.5626E-01 -1.9867E-01 -5.4093E-01  3.8528E-01  1.6760E-01 -7.9309E-02 -5.0187E+00  7.2590E-01  2.3812E-01
             1.8843E-01
 GRADIENT:   3.9092E+02  6.6263E+02  5.5315E-01  9.8974E+01  1.6358E+01  8.9259E+01  5.7895E+00  0.0000E+00  3.5714E+01  1.4469E+00
             1.6094E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1630.10510120802        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1636
 NPARAMETR:  9.9646E-01  1.7435E+00  7.4254E-01  5.2724E-01  1.3297E+00  1.0699E+00  8.3611E-01  1.0000E-02  1.8684E+00  1.1479E+00
             1.0924E+00
 PARAMETER:  9.6457E-02  6.5591E-01 -1.9767E-01 -5.4010E-01  3.8493E-01  1.6760E-01 -7.8998E-02 -5.0187E+00  7.2510E-01  2.3794E-01
             1.8837E-01
 GRADIENT:   1.6113E+00 -9.3292E+00 -1.9516E-01 -3.6951E-01  3.4250E-01  5.8210E-01 -2.8837E-02  0.0000E+00  3.4429E-01  3.3315E-02
             5.4397E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1630.10605615946        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1830             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9647E-01  1.7424E+00  7.4403E-01  5.2807E-01  1.3287E+00  1.0699E+00  8.3665E-01  1.0000E-02  1.8655E+00  1.1474E+00
             1.0922E+00
 PARAMETER:  9.6459E-02  6.5525E-01 -1.9568E-01 -5.3852E-01  3.8417E-01  1.6761E-01 -7.8349E-02 -5.0187E+00  7.2350E-01  2.3749E-01
             1.8822E-01
 GRADIENT:   3.9111E+02  6.6122E+02  1.0467E+00  9.8854E+01  1.5411E+01  8.9298E+01  5.8332E+00  0.0000E+00  3.5437E+01  1.4002E+00
             1.4433E+00

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1630.10605615946        NO. OF FUNC. EVALS.:  69
 CUMULATIVE NO. OF FUNC. EVALS.:     1899
 NPARAMETR:  9.9645E-01  1.7406E+00  7.4257E-01  5.2908E-01  1.3293E+00  1.0699E+00  8.3673E-01  1.0000E-02  1.8659E+00  1.1474E+00
             1.0926E+00
 PARAMETER:  9.6459E-02  6.5525E-01 -1.9568E-01 -5.3852E-01  3.8417E-01  1.6761E-01 -7.8349E-02 -5.0187E+00  7.2350E-01  2.3749E-01
             1.8822E-01
 GRADIENT:   1.2397E-02  8.4984E-01  1.5435E-01 -3.1688E-01 -2.1010E-01  2.6799E-03 -5.4932E-03  0.0000E+00 -1.9469E-02 -1.3738E-03
            -6.1014E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1899
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.4724E-05 -3.8870E-02 -2.3746E-04  2.8427E-02 -4.4569E-02
 SE:             2.9828E-02  2.1634E-02  9.0653E-05  2.3198E-02  2.2432E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9854E-01  7.2378E-02  8.8080E-03  2.2041E-01  4.6942E-02

 ETASHRINKSD(%)  7.2801E-02  2.7524E+01  9.9696E+01  2.2285E+01  2.4849E+01
 ETASHRINKVR(%)  1.4555E-01  4.7472E+01  9.9999E+01  3.9604E+01  4.3524E+01
 EBVSHRINKSD(%)  4.6257E-01  2.5442E+01  9.9723E+01  2.4339E+01  2.3088E+01
 EBVSHRINKVR(%)  9.2300E-01  4.4411E+01  9.9999E+01  4.2754E+01  4.0846E+01
 RELATIVEINF(%)  9.9010E+01  5.4702E+00  1.9521E-04  6.4953E+00  2.0514E+01
 EPSSHRINKSD(%)  4.2846E+01
 EPSSHRINKVR(%)  6.7334E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1630.1060561594622     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -894.95522959572406     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.16
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.06
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1630.106       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.96E-01  1.74E+00  7.44E-01  5.28E-01  1.33E+00  1.07E+00  8.37E-01  1.00E-02  1.87E+00  1.15E+00  1.09E+00
 


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
+        9.69E+02
 
 TH 2
+       -5.19E+00  3.15E+02
 
 TH 3
+        4.97E+00  9.44E+01  1.61E+02
 
 TH 4
+       -9.16E+00  2.72E+02 -1.10E+02  6.76E+02
 
 TH 5
+       -2.74E+00 -1.08E+02 -1.19E+02  1.06E+02  2.85E+02
 
 TH 6
+        4.08E-01 -6.13E-01  1.14E+00 -2.94E+00 -7.03E-01  1.71E+02
 
 TH 7
+        9.40E-01  5.33E+00  1.35E+01 -2.21E+01 -1.68E+01 -3.22E-01  9.79E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.20E-01 -1.77E+01 -1.83E+01  4.39E+01  1.60E+00 -3.81E-01  1.39E+01  0.00E+00  2.56E+01
 
 TH10
+        2.88E-01 -9.30E+00 -1.84E+01 -3.52E+00 -4.53E+01 -9.19E-03  7.69E+00  0.00E+00  3.12E+00  5.93E+01
 
 TH11
+       -7.88E+00 -2.17E+01 -2.68E+01  3.17E+00  1.51E-01  2.80E+00  1.08E+01  0.00E+00  2.47E+00  1.92E+01  1.83E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       36.270
Stop Time:
Wed Sep 29 18:59:19 CDT 2021
