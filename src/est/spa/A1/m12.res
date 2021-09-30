Wed Sep 29 11:55:25 CDT 2021
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
$DATA ../../../../data/spa/A1/dat12.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m12.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1388.16730593263        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3090E+02 -1.3615E+01  6.7569E+00 -1.0330E+01  4.2453E+01  5.2835E+01 -2.8899E+01  8.9137E+00 -2.2859E+01 -2.8274E+01
            -4.5821E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1486.30807895093        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0967E+00  1.0823E+00  1.0508E+00  1.0781E+00  1.0214E+00  1.1088E+00  1.1622E+00  8.8523E-01  1.0394E+00  1.0357E+00
             2.0276E+00
 PARAMETER:  1.9233E-01  1.7904E-01  1.4958E-01  1.7519E-01  1.2115E-01  2.0330E-01  2.5028E-01 -2.1911E-02  1.3869E-01  1.3512E-01
             8.0685E-01
 GRADIENT:   3.5421E+02  5.4217E+01 -1.5072E+00  8.7945E+01 -1.2018E+01  5.2238E+01 -1.1360E+00  6.3746E+00  9.1299E+00  2.2041E+00
             6.2638E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1500.24336624263        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0435E+00  8.6013E-01  8.4708E-01  1.1713E+00  8.6913E-01  1.0925E+00  1.1818E+00  9.4950E-02  9.3586E-01  1.0323E+00
             1.7257E+00
 PARAMETER:  1.4259E-01 -5.0674E-02 -6.5962E-02  2.5815E-01 -4.0264E-02  1.8847E-01  2.6701E-01 -2.2544E+00  3.3705E-02  1.3180E-01
             6.4563E-01
 GRADIENT:   2.7237E+02  6.1021E+00 -4.2703E+01  1.2682E+02  5.2086E+01  6.3493E+01 -1.6784E+01  1.8396E-01 -3.2203E+00  7.8272E+00
             4.0894E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1507.64817664962        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      250
 NPARAMETR:  9.8254E-01  8.5405E-01  7.0878E-01  1.1143E+00  7.6038E-01  9.7084E-01  1.4579E+00  1.1631E-01  9.4410E-01  8.6873E-01
             1.6224E+00
 PARAMETER:  8.2391E-02 -5.7765E-02 -2.4421E-01  2.0821E-01 -1.7393E-01  7.0412E-02  4.7698E-01 -2.0515E+00  4.2473E-02 -4.0724E-02
             5.8392E-01
 GRADIENT:  -3.3092E+01 -4.3725E+00 -2.8780E+01 -1.7018E+00  2.3656E+01 -3.3415E-01 -2.8188E+00  3.5868E-01  9.0638E+00  7.8056E+00
            -2.3744E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1515.95882306702        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      428
 NPARAMETR:  9.8969E-01  4.5485E-01  8.2779E-01  1.3746E+00  6.8138E-01  9.5350E-01  2.4602E+00  5.9191E-02  7.8192E-01  8.2061E-01
             1.7387E+00
 PARAMETER:  8.9640E-02 -6.8778E-01 -8.8993E-02  4.1820E-01 -2.8364E-01  5.2388E-02  1.0002E+00 -2.7270E+00 -1.4600E-01 -9.7711E-02
             6.5314E-01
 GRADIENT:  -6.8582E+00  1.4453E+01  2.4643E+00  2.7581E+01 -4.3103E+00 -3.7535E+00  4.3159E+00  6.4277E-02 -5.0057E+00 -8.4456E-01
             6.3964E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1518.56555354453        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      603
 NPARAMETR:  9.8926E-01  2.2073E-01  8.3886E-01  1.4810E+00  6.3134E-01  9.6430E-01  3.3237E+00  1.0000E-02  7.9099E-01  8.5147E-01
             1.7074E+00
 PARAMETER:  8.9202E-02 -1.4108E+00 -7.5713E-02  4.9274E-01 -3.5992E-01  6.3644E-02  1.3011E+00 -4.9786E+00 -1.3447E-01 -6.0794E-02
             6.3496E-01
 GRADIENT:   5.3068E+00  7.6166E-01  5.9326E+00 -4.8676E+00 -7.3969E+00  1.3556E+00 -1.2787E+00  0.0000E+00  1.2911E+00 -1.2813E-01
             4.3912E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1518.85731226093        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      779
 NPARAMETR:  9.8349E-01  1.2666E-01  8.6375E-01  1.5473E+00  6.2516E-01  9.5996E-01  4.0984E+00  1.0000E-02  7.9064E-01  8.9136E-01
             1.7024E+00
 PARAMETER:  8.3350E-02 -1.9663E+00 -4.6471E-02  5.3650E-01 -3.6974E-01  5.9134E-02  1.5106E+00 -6.7782E+00 -1.3491E-01 -1.5010E-02
             6.3205E-01
 GRADIENT:  -2.7857E+00  1.5880E-01  5.8135E+00  2.0696E+01 -1.0447E+01 -2.2679E-01 -4.3745E+00  0.0000E+00  3.6218E+00  1.8328E+00
            -5.2047E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1519.20661994601        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      956
 NPARAMETR:  9.8360E-01  8.2090E-02  8.8458E-01  1.5678E+00  6.2921E-01  9.6104E-01  4.9327E+00  1.0000E-02  7.7328E-01  9.0564E-01
             1.7049E+00
 PARAMETER:  8.3460E-02 -2.3999E+00 -2.2644E-02  5.4965E-01 -3.6328E-01  6.0263E-02  1.6959E+00 -8.1623E+00 -1.5711E-01  8.8677E-04
             6.3350E-01
 GRADIENT:   7.3805E-01  5.9166E-01  1.5859E+00  3.8764E+00 -1.3317E+00  4.8341E-01  1.3680E+00  0.0000E+00  8.5594E-02 -1.1351E+00
            -6.6102E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1519.23048169336        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1131
 NPARAMETR:  9.8331E-01  8.1678E-02  8.7054E-01  1.5627E+00  6.2279E-01  9.5985E-01  4.9013E+00  1.0000E-02  7.7230E-01  9.0327E-01
             1.7043E+00
 PARAMETER:  8.3170E-02 -2.4050E+00 -3.8646E-02  5.4639E-01 -3.7354E-01  5.9024E-02  1.6895E+00 -8.2315E+00 -1.5839E-01 -1.7386E-03
             6.3314E-01
 GRADIENT:   1.8683E-02  9.6261E-02 -1.1371E-01 -3.5085E-01  1.5555E-01 -9.7354E-03  1.7444E-01  0.0000E+00 -1.6089E-01 -9.2329E-03
            -3.6169E-03

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1519.23048169336        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1153
 NPARAMETR:  9.8331E-01  8.1678E-02  8.7054E-01  1.5627E+00  6.2279E-01  9.5985E-01  4.9013E+00  1.0000E-02  7.7230E-01  9.0327E-01
             1.7043E+00
 PARAMETER:  8.3170E-02 -2.4050E+00 -3.8646E-02  5.4639E-01 -3.7354E-01  5.9024E-02  1.6895E+00 -8.2315E+00 -1.5839E-01 -1.7386E-03
             6.3314E-01
 GRADIENT:   1.8683E-02  9.6261E-02 -1.1371E-01 -3.5085E-01  1.5555E-01 -9.7354E-03  1.7444E-01  0.0000E+00 -1.6089E-01 -9.2329E-03
            -3.6169E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1153
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.6779E-04  3.1862E-02 -1.7013E-04 -2.6796E-02 -8.2435E-03
 SE:             2.9500E-02  1.4389E-02  1.6807E-04  2.6363E-02  2.1639E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7653E-01  2.6811E-02  3.1142E-01  3.0943E-01  7.0323E-01

 ETASHRINKSD(%)  1.1706E+00  5.1794E+01  9.9437E+01  1.1680E+01  2.7508E+01
 ETASHRINKVR(%)  2.3274E+00  7.6761E+01  9.9997E+01  2.1995E+01  4.7449E+01
 EBVSHRINKSD(%)  1.2370E+00  6.9338E+01  9.9354E+01  8.4944E+00  2.1453E+01
 EBVSHRINKVR(%)  2.4586E+00  9.0598E+01  9.9996E+01  1.6267E+01  3.8303E+01
 RELATIVEINF(%)  9.7174E+01  4.5506E+00  3.0206E-04  3.9758E+01  4.4460E+00
 EPSSHRINKSD(%)  3.7933E+01
 EPSSHRINKVR(%)  6.1477E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1519.2304816933579     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -784.07965512961971     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.44
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.15
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1519.230       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.83E-01  8.17E-02  8.71E-01  1.56E+00  6.23E-01  9.60E-01  4.90E+00  1.00E-02  7.72E-01  9.03E-01  1.70E+00
 


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
+        1.22E+03
 
 TH 2
+        2.98E+01  3.00E+04
 
 TH 3
+       -1.38E+01 -7.88E+02  8.59E+02
 
 TH 4
+       -3.90E+01 -8.40E+02  1.21E+02  8.95E+02
 
 TH 5
+        4.78E+01  1.47E+03 -1.56E+03 -4.86E+02  3.13E+03
 
 TH 6
+       -2.07E-03 -4.22E+01  1.25E+01  1.56E+00 -1.72E+01  2.06E+02
 
 TH 7
+        2.22E+00  8.14E+02 -4.40E+01 -1.86E+02  5.93E+02 -1.28E+00  2.51E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -3.17E+00 -7.65E+02  1.45E+02  1.09E+02 -2.54E+02  3.88E+00 -1.03E+03  0.00E+00  3.45E+02
 
 TH10
+       -1.70E+01 -6.75E+02  1.81E+02  1.66E+02 -4.65E+02  6.33E+00 -1.53E+01  0.00E+00  1.00E+02  2.67E+02
 
 TH11
+       -1.35E+01 -7.80E+01 -2.84E+00  4.60E+00 -2.67E+01  3.94E+00 -1.22E+00  0.00E+00  2.33E+01  3.76E+01  8.79E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.655
Stop Time:
Wed Sep 29 11:55:48 CDT 2021
