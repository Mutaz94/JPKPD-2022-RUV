Sat Sep 25 01:20:23 CDT 2021
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
$DATA ../../../../data/int/SL2/dat53.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m53.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2038.43166606026        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.1515E+01 -1.1236E+01  1.0890E+02 -1.0400E+01  1.1324E+02 -3.0315E+01 -8.0754E+01 -1.8664E+02 -4.1432E+01 -1.9121E+01
            -3.3161E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3026.70178702897        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.5751E-01  1.3488E+00  1.1254E+00  8.5491E-01  1.2309E+00  1.0795E+00  8.8959E-01  8.8289E-01  7.9813E-01  1.0998E+00
             1.9718E+00
 PARAMETER:  5.6582E-02  3.9923E-01  2.1814E-01 -5.6762E-02  3.0777E-01  1.7647E-01 -1.6995E-02 -2.4559E-02 -1.2548E-01  1.9516E-01
             7.7893E-01
 GRADIENT:  -4.2875E+01  1.1763E+01  8.6109E+00  7.6667E-01 -2.2152E+00  1.6393E+00  3.9803E-01 -2.4764E+00 -2.7792E+01 -2.6189E+01
            -1.8331E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3036.52646909518        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.8495E-01  1.5601E+00  1.1085E+00  7.6202E-01  1.3500E+00  1.1223E+00  8.1397E-01  3.0924E-01  1.0266E+00  1.2127E+00
             2.0525E+00
 PARAMETER:  8.4840E-02  5.4473E-01  2.0297E-01 -1.7179E-01  4.0013E-01  2.1536E-01 -1.0583E-01 -1.0736E+00  1.2629E-01  2.9285E-01
             8.1907E-01
 GRADIENT:   7.8972E+00  8.9160E+01  1.0921E+01  4.5397E+01 -1.9881E+01  1.7832E+01  5.9480E+00 -4.0001E-01 -2.7540E+00 -1.8315E+01
            -8.4920E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3041.33263702682        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.8137E-01  1.5927E+00  1.0434E+00  6.9258E-01  1.4081E+00  1.0644E+00  6.9684E-01  1.6907E-01  1.1061E+00  1.3436E+00
             2.1113E+00
 PARAMETER:  8.1194E-02  5.6545E-01  1.4249E-01 -2.6734E-01  4.4228E-01  1.6237E-01 -2.6119E-01 -1.6774E+00  2.0085E-01  3.9536E-01
             8.4731E-01
 GRADIENT:  -4.4429E-01 -1.0327E+01  1.0490E+00 -4.5160E+00 -3.1302E+00 -2.5797E+00 -4.3076E+00 -5.4412E-02  8.6615E-01  2.0421E+00
            -1.1732E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3041.99146192660        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  9.8161E-01  1.6291E+00  1.0338E+00  6.7702E-01  1.4375E+00  1.0712E+00  7.3032E-01  1.6030E-01  1.0877E+00  1.3430E+00
             2.1222E+00
 PARAMETER:  8.1439E-02  5.8801E-01  1.3327E-01 -2.9005E-01  4.6290E-01  1.6879E-01 -2.1427E-01 -1.7307E+00  1.8408E-01  3.9493E-01
             8.5246E-01
 GRADIENT:  -3.6016E-01  6.3482E+00  3.5553E-01  3.0249E+00  6.9916E-01 -9.1099E-02  1.4103E+00 -4.5485E-02 -1.2337E+00 -4.7679E-01
             1.3650E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3042.02575958114        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  9.8174E-01  1.6319E+00  1.0257E+00  6.7025E-01  1.4384E+00  1.0712E+00  7.1151E-01  1.8177E-01  1.1142E+00  1.3480E+00
             2.1187E+00
 PARAMETER:  8.1567E-02  5.8975E-01  1.2539E-01 -3.0011E-01  4.6351E-01  1.6879E-01 -2.4037E-01 -1.6050E+00  2.0813E-01  3.9860E-01
             8.5080E-01
 GRADIENT:   1.2487E-01 -3.2822E+00 -6.8607E-02 -1.7323E+00 -3.9158E-01 -3.6020E-02 -8.1738E-01 -5.4250E-02  3.5835E-01  1.6712E-01
            -1.9684E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3044.78253744623        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      555
 NPARAMETR:  9.8661E-01  2.0450E+00  8.5916E-01  4.3089E-01  1.7559E+00  1.0756E+00  6.2574E-01  8.2357E-01  1.4932E+00  1.5250E+00
             2.1060E+00
 PARAMETER:  8.6523E-02  8.1542E-01 -5.1802E-02 -7.4190E-01  6.6298E-01  1.7287E-01 -3.6881E-01 -9.4105E-02  5.0092E-01  5.2200E-01
             8.4481E-01
 GRADIENT:  -5.8398E-01  4.9106E+01  3.6437E-01  1.9433E+01  4.1303E+00 -1.7843E+00 -3.9919E+00 -1.3345E-03 -1.4202E+00  1.9347E+00
            -8.8452E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3046.21278479248        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      734
 NPARAMETR:  9.8595E-01  2.2467E+00  6.3301E-01  2.8580E-01  1.8807E+00  1.0789E+00  5.9904E-01  3.8199E-01  1.9953E+00  1.5781E+00
             2.1166E+00
 PARAMETER:  8.5848E-02  9.0946E-01 -3.5727E-01 -1.1525E+00  7.3165E-01  1.7597E-01 -4.1243E-01 -8.6235E-01  7.9077E-01  5.5622E-01
             8.4982E-01
 GRADIENT:  -1.6977E+00  3.4012E+01  4.0535E-01  8.5085E+00 -2.4026E+00 -5.3232E-01 -2.8161E+00  1.1361E-01  1.1528E+00 -8.2057E-01
             6.9449E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3046.90006177696        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      912
 NPARAMETR:  9.8592E-01  2.3587E+00  4.8696E-01  1.9861E-01  1.9819E+00  1.0830E+00  5.9052E-01  5.7558E-02  2.4408E+00  1.6290E+00
             2.1082E+00
 PARAMETER:  8.5818E-02  9.5810E-01 -6.1957E-01 -1.5164E+00  7.8407E-01  1.7976E-01 -4.2675E-01 -2.7550E+00  9.9233E-01  5.8800E-01
             8.4583E-01
 GRADIENT:  -1.2122E+00  1.3143E+00 -1.7281E+00  2.7161E+00  3.2403E+00  9.8532E-01 -1.3243E+00  4.7656E-03  1.0822E+00 -2.8074E-01
             3.4981E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3047.15790058232        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1088
 NPARAMETR:  9.8686E-01  2.4601E+00  4.3389E-01  1.3131E-01  2.0621E+00  1.0806E+00  5.8258E-01  1.0000E-02  3.0634E+00  1.6732E+00
             2.1056E+00
 PARAMETER:  8.6774E-02  1.0002E+00 -7.3496E-01 -1.9302E+00  8.2373E-01  1.7751E-01 -4.4028E-01 -4.8582E+00  1.2195E+00  6.1476E-01
             8.4460E-01
 GRADIENT:   6.2687E-01  1.1933E+01 -6.9649E-01  1.0716E+00  1.9602E+00  1.1881E-01 -1.0995E+00  0.0000E+00 -7.0474E-01  6.5546E-02
            -1.1900E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3047.19534710414        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1263
 NPARAMETR:  9.8650E-01  2.4561E+00  4.4668E-01  1.2965E-01  2.0579E+00  1.0801E+00  5.8504E-01  1.0000E-02  3.1143E+00  1.6715E+00
             2.1058E+00
 PARAMETER:  8.6404E-02  9.9856E-01 -7.0591E-01 -1.9429E+00  8.2171E-01  1.7704E-01 -4.3608E-01 -4.9772E+00  1.2360E+00  6.1369E-01
             8.4468E-01
 GRADIENT:   5.8107E-02  1.3617E-01 -2.8935E-03  1.5709E-02 -4.0184E-02 -2.6417E-02  3.9696E-03  0.0000E+00 -2.1149E-04 -4.7629E-03
            -4.6710E-03

0ITERATION NO.:   52    OBJECTIVE VALUE:  -3047.19534935427        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:     1326
 NPARAMETR:  9.8642E-01  2.4559E+00  4.4661E-01  1.2948E-01  2.0583E+00  1.0803E+00  5.8501E-01  1.0000E-02  3.1161E+00  1.6716E+00
             2.1058E+00
 PARAMETER:  8.6388E-02  9.9862E-01 -7.0627E-01 -1.9439E+00  8.2180E-01  1.7708E-01 -4.3610E-01 -4.9836E+00  1.2365E+00  6.1374E-01
             8.4468E-01
 GRADIENT:   2.6249E-02  1.5960E-01 -1.4375E-03  6.2789E-03 -1.6556E-02 -1.5179E-02  2.0617E-03  0.0000E+00 -1.7613E-03 -3.5418E-03
            -3.8446E-05

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1326
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0546E-03 -3.0326E-02 -2.4103E-05  4.0851E-02 -1.9181E-02
 SE:             2.9652E-02  2.5374E-02  1.9667E-05  1.6772E-02  2.7008E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7163E-01  2.3203E-01  2.2037E-01  1.4864E-02  4.7757E-01

 ETASHRINKSD(%)  6.6221E-01  1.4993E+01  9.9934E+01  4.3813E+01  9.5211E+00
 ETASHRINKVR(%)  1.3200E+00  2.7738E+01  1.0000E+02  6.8430E+01  1.8136E+01
 EBVSHRINKSD(%)  9.1022E-01  1.2869E+01  9.9908E+01  5.2942E+01  6.3015E+00
 EBVSHRINKVR(%)  1.8122E+00  2.4082E+01  1.0000E+02  7.7855E+01  1.2206E+01
 RELATIVEINF(%)  9.8164E+01  1.9115E+01  5.8906E-05  5.3298E+00  6.4431E+01
 EPSSHRINKSD(%)  1.6843E+01
 EPSSHRINKVR(%)  3.0850E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3047.1953493542742     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1393.1059895858634     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.38
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.20
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3047.195       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.86E-01  2.46E+00  4.47E-01  1.30E-01  2.06E+00  1.08E+00  5.85E-01  1.00E-02  3.12E+00  1.67E+00  2.11E+00
 


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
+        9.58E+02
 
 TH 2
+       -9.74E+00  4.47E+02
 
 TH 3
+       -2.12E+00  4.58E+01  9.00E+01
 
 TH 4
+       -6.16E+00  4.10E+02 -2.27E+02  2.36E+03
 
 TH 5
+       -1.71E+00 -1.74E+01 -5.47E+00  5.13E+01  8.33E+01
 
 TH 6
+        6.48E+00 -3.39E+00 -1.46E-01 -3.52E+00 -8.58E-01  1.64E+02
 
 TH 7
+        1.10E+00 -6.89E+00 -1.16E+01 -2.95E+01 -3.08E+00 -1.58E+00  4.08E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.90E-01 -1.11E+01 -1.14E+01  9.58E+01 -1.16E+00  8.87E-02 -6.36E-01  0.00E+00  7.39E+00
 
 TH10
+       -1.01E+00 -4.69E+00  6.06E+00  5.75E+00 -4.98E+00 -1.00E-01 -5.70E+00  0.00E+00 -2.50E-01  5.40E+01
 
 TH11
+       -1.06E+01 -1.24E+01  8.56E+00 -4.13E+01  1.41E+00  1.66E+00  1.00E+01  0.00E+00 -3.45E-01  5.13E+00  2.68E+02
 
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
 #CPUT: Total CPU Time in Seconds,       43.685
Stop Time:
Sat Sep 25 01:21:08 CDT 2021
