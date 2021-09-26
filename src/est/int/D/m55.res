Sat Sep 25 05:58:31 CDT 2021
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
$DATA ../../../../data/int/D/dat55.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m55.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   37221.8420354251        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7160E+02  5.5520E+02 -3.5498E+01  3.6492E+02  1.6475E+02 -2.5469E+03 -1.1068E+03 -9.5535E+01 -1.8092E+03 -8.0271E+02
            -7.5357E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -755.774439516911        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.2021E+00  1.8407E+00  8.9293E-01  2.3488E+00  9.5139E-01  4.9100E+00  3.3447E+00  9.7919E-01  2.0220E+00  1.4750E+00
             1.3098E+01
 PARAMETER:  2.8407E-01  7.1014E-01 -1.3252E-02  9.5391E-01  5.0171E-02  1.6913E+00  1.3074E+00  7.8969E-02  8.0408E-01  4.8867E-01
             2.6724E+00
 GRADIENT:  -1.3904E+01  3.3611E+01 -4.2826E+01  1.9421E+02  6.7244E+00  1.3898E+02 -4.2435E+01  4.9234E+00 -4.0062E+01  3.4816E+01
             2.7041E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -815.351691486905        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0485E+00  1.6523E+00  3.0055E+01  6.1535E+00  2.8041E+00  2.5397E+00  8.0259E+00  6.8522E-01  4.5971E+00  1.5510E+00
             1.3054E+01
 PARAMETER:  1.4731E-01  6.0216E-01  3.5030E+00  1.9170E+00  1.1311E+00  1.0321E+00  2.1827E+00 -2.7801E-01  1.6254E+00  5.3893E-01
             2.6691E+00
 GRADIENT:  -4.5322E+01  2.1449E+01 -1.2624E+01  1.0127E+02  2.9084E+01  7.3037E+01  1.2801E+01  4.0938E-01 -5.6632E-01  2.9206E+01
             3.0830E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -955.243541671046        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.2591E+00  1.1869E+00  1.0168E+01  1.7374E+00  1.8593E+00  2.7735E+00  5.3090E+00  1.5427E+00  1.3888E+00  7.0082E-01
             1.2618E+01
 PARAMETER:  3.3039E-01  2.7133E-01  2.4192E+00  6.5242E-01  7.2022E-01  1.1201E+00  1.7694E+00  5.3351E-01  4.2846E-01 -2.5550E-01
             2.6351E+00
 GRADIENT:  -1.2050E+01  1.3964E+01  2.5216E+00  3.6479E+01 -3.5418E+01  5.7756E+01 -7.8815E+00  6.6326E-01  5.4651E+00  6.8721E+00
             3.4104E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1025.23101446305        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.0909E+00  7.0907E-01  2.1784E+01  1.5300E+00  2.3184E+00  1.9169E+00  5.8668E+00  1.4185E+00  1.2675E+00  5.4522E-01
             1.0074E+01
 PARAMETER:  1.8698E-01 -2.4380E-01  3.1812E+00  5.2530E-01  9.4088E-01  7.5068E-01  1.8693E+00  4.4957E-01  3.3708E-01 -5.0656E-01
             2.4100E+00
 GRADIENT:  -1.0697E+01 -5.7754E+00 -1.4788E+00 -2.4035E+01  7.4600E+00  8.2082E+00  3.2515E+00  7.1700E-02  7.9934E+00  3.3425E+00
             6.0119E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1028.09753722893        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.0921E+00  6.3078E-01  3.0384E+01  1.6080E+00  2.3575E+00  1.8514E+00  6.0790E+00  1.4843E+00  1.1949E+00  2.7700E-01
             9.8202E+00
 PARAMETER:  1.8810E-01 -3.6079E-01  3.5139E+00  5.7498E-01  9.5760E-01  7.1597E-01  1.9048E+00  4.9492E-01  2.7807E-01 -1.1837E+00
             2.3844E+00
 GRADIENT:  -4.8199E+00 -2.7959E+00 -4.6380E-01  1.4022E-01 -1.0565E+00 -1.2310E+00 -1.1735E+00  1.9851E-02 -1.0242E+00  6.7559E-01
            -2.1774E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1028.45449920534        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      468
 NPARAMETR:  1.0983E+00  7.2740E-01  4.2405E+01  1.5656E+00  2.4337E+00  1.8582E+00  5.8801E+00  1.7853E+00  1.1671E+00  9.8021E-02
             9.8731E+00
 PARAMETER:  1.9377E-01 -2.1828E-01  3.8473E+00  5.4824E-01  9.8942E-01  7.1960E-01  1.8716E+00  6.7961E-01  2.5455E-01 -2.2226E+00
             2.3898E+00
 GRADIENT:  -3.6110E+00 -1.6122E+00 -2.8815E-01 -1.7309E+00  8.4736E-01 -1.0750E+00 -1.0662E+01  1.2478E-02 -4.9062E-01  8.4280E-02
            -2.9356E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1029.25246542839        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      646
 NPARAMETR:  1.1090E+00  5.6431E-01  6.5808E+01  1.6906E+00  2.4894E+00  1.8584E+00  6.7581E+00  2.0221E+00  1.2451E+00  8.5148E-02
             9.8550E+00
 PARAMETER:  2.0350E-01 -4.7216E-01  4.2867E+00  6.2510E-01  1.0120E+00  7.1974E-01  2.0107E+00  8.0414E-01  3.1923E-01 -2.3634E+00
             2.3880E+00
 GRADIENT:   2.5242E+00  3.6669E-01 -2.5098E-01 -1.6595E+00  3.7398E+00 -4.7453E-01  4.2397E-01 -9.5039E-04  3.8974E-02  6.2578E-02
            -4.9063E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1029.38749620970        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      821
 NPARAMETR:  1.1026E+00  5.0832E-01  1.4399E+02  1.7336E+00  2.5156E+00  1.8532E+00  6.9638E+00  2.8005E+00  1.2760E+00  1.7883E-02
             9.8770E+00
 PARAMETER:  1.9765E-01 -5.7665E-01  5.0697E+00  6.5018E-01  1.0225E+00  7.1694E-01  2.0407E+00  1.1298E+00  3.4377E-01 -3.9239E+00
             2.3902E+00
 GRADIENT:  -1.7447E-01 -7.3388E-02 -1.3133E-02 -2.6050E-01  9.9385E-02  1.0029E-02 -8.6721E-02 -2.0215E-03  8.6044E-02  2.6847E-03
             4.2615E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1029.39052213750        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1002             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1027E+00  5.0944E-01  2.0459E+02  1.7346E+00  2.5279E+00  1.8523E+00  6.9709E+00  3.2792E+00  1.2752E+00  1.0000E-02
             9.8766E+00
 PARAMETER:  1.9776E-01 -5.7445E-01  5.4210E+00  6.5077E-01  1.0274E+00  7.1643E-01  2.0417E+00  1.2876E+00  3.4309E-01 -4.7480E+00
             2.3902E+00
 GRADIENT:   9.8932E-01  3.9372E-01  1.1869E-03  2.2121E+00  3.2907E-01  1.2434E+00  1.3484E+01 -1.3672E-03  1.5376E-01  0.0000E+00
             4.2852E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1029.39054617192        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1164
 NPARAMETR:  1.1027E+00  5.0881E-01  2.0119E+02  1.7348E+00  2.5277E+00  1.8522E+00  6.9705E+00  3.2845E+00  1.2751E+00  1.0000E-02
             9.8767E+00
 PARAMETER:  1.9776E-01 -5.7569E-01  5.4042E+00  6.5090E-01  1.0273E+00  7.1638E-01  2.0417E+00  1.2892E+00  3.4300E-01 -4.7480E+00
             2.3902E+00
 GRADIENT:   9.7101E-03 -6.2934E-03 -2.0061E-04 -3.5872E-02  5.3678E-02  3.6808E-03  3.9102E-02 -1.4728E-03 -1.5820E-02  0.0000E+00
             2.9220E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1164
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7164E-02  5.6411E-02  3.5976E-05 -8.5199E-02 -2.4353E-06
 SE:             2.8651E-02  2.1843E-02  3.8057E-04  1.6097E-02  1.2263E-04
 N:                     100         100         100         100         100

 P VAL.:         5.4913E-01  9.8065E-03  9.2469E-01  1.2073E-07  9.8416E-01

 ETASHRINKSD(%)  4.0148E+00  2.6824E+01  9.8725E+01  4.6072E+01  9.9589E+01
 ETASHRINKVR(%)  7.8683E+00  4.6453E+01  9.9984E+01  7.0918E+01  9.9998E+01
 EBVSHRINKSD(%)  6.6937E+00  2.5451E+01  9.8444E+01  4.0863E+01  9.9549E+01
 EBVSHRINKVR(%)  1.2939E+01  4.4425E+01  9.9976E+01  6.5028E+01  9.9998E+01
 RELATIVEINF(%)  8.6831E+01  2.8053E+01  4.2911E-03  1.7357E+01  3.5201E-04
 EPSSHRINKSD(%)  5.6361E+00
 EPSSHRINKVR(%)  1.0955E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1029.3905461719150     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       624.69881359649571     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    35.01
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.92
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1029.391       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.10E+00  5.09E-01  2.01E+02  1.73E+00  2.53E+00  1.85E+00  6.97E+00  3.28E+00  1.28E+00  1.00E-02  9.88E+00
 


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
+        1.93E+02
 
 TH 2
+        6.19E+00  4.79E+01
 
 TH 3
+       -4.89E-03  1.35E-03  1.06E-06
 
 TH 4
+       -2.73E+00  3.75E+01  8.47E-05  1.01E+02
 
 TH 5
+       -7.63E-02 -7.10E+00 -4.53E-03 -1.15E+01  3.19E+01
 
 TH 6
+       -1.81E+01  1.17E-01 -4.92E-04  6.27E-01  1.89E+00  2.86E+01
 
 TH 7
+        7.82E-01  5.37E+00 -4.52E-05 -4.14E+00  2.04E-01  4.22E-01  1.92E+00
 
 TH 8
+        1.77E-01  8.11E-01  2.14E-04 -5.00E-02 -2.01E-02  1.75E-01  4.21E-03  2.38E-02
 
 TH 9
+        2.00E+00  1.64E+00  1.55E-03 -3.19E+01  3.58E+00  7.30E-01  1.36E+00 -4.19E-01  2.74E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.27E+00 -4.55E+00 -6.11E-05 -6.98E+00  2.76E-01  3.55E+00  1.01E-03 -4.28E-03  2.95E+00  0.00E+00  9.80E+00
 
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
 #CPUT: Total CPU Time in Seconds,       53.064
Stop Time:
Sat Sep 25 05:59:26 CDT 2021
