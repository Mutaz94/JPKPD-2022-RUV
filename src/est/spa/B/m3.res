Sat Sep 18 08:15:26 CDT 2021
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
$DATA ../../../../data/spa/B/dat3.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m3.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1718.52591078141        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -4.4718E+00 -7.8752E+01 -5.0790E+01 -5.2188E+01  4.4189E+01  5.6462E+00  4.1961E+00  1.1413E+01  1.2610E+01  3.0494E+01
             2.5419E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1727.08945162642        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  1.0029E+00  1.1013E+00  1.1522E+00  9.9902E-01  1.0459E+00  9.8182E-01  9.3940E-01  9.4543E-01  9.4741E-01  7.7984E-01
             1.0353E+00
 PARAMETER:  1.0293E-01  1.9651E-01  2.4167E-01  9.9021E-02  1.4492E-01  8.1657E-02  3.7487E-02  4.3882E-02  4.5977E-02 -1.4866E-01
             1.3466E-01
 GRADIENT:  -1.8820E+00  1.8281E+01  1.2303E+01  1.4779E+01  8.0280E+00 -1.8274E+00 -7.3185E-01 -6.3689E+00 -6.8065E+00 -7.6607E+00
             2.4145E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1728.14068824376        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      250
 NPARAMETR:  1.0262E+00  1.2037E+00  1.0427E+00  9.3765E-01  1.0372E+00  9.9746E-01  9.3237E-01  1.0013E+00  9.9761E-01  7.1961E-01
             1.0271E+00
 PARAMETER:  1.2585E-01  2.8538E-01  1.4178E-01  3.5617E-02  1.3649E-01  9.7455E-02  2.9979E-02  1.0130E-01  9.7602E-02 -2.2904E-01
             1.2675E-01
 GRADIENT:   4.3607E+00  2.2093E+01  1.1783E+01  1.6862E+01 -3.5043E-01 -2.0285E+00  1.6611E+00 -3.6319E+00 -3.8838E+00 -9.1008E+00
            -6.8493E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1729.64476110614        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      426
 NPARAMETR:  1.0256E+00  1.4117E+00  7.0855E-01  7.8687E-01  9.8907E-01  1.0124E+00  8.4825E-01  6.5363E-01  1.1445E+00  7.7562E-01
             1.0037E+00
 PARAMETER:  1.2532E-01  4.4481E-01 -2.4453E-01 -1.3969E-01  8.9008E-02  1.1231E-01 -6.4583E-02 -3.2521E-01  2.3495E-01 -1.5409E-01
             1.0365E-01
 GRADIENT:  -1.1590E+00  1.5458E+01 -1.2197E+00  1.4536E+01 -2.2350E+01  2.5416E+00  4.6605E+00  1.8353E+00  4.4143E+00  8.2749E+00
            -1.3489E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1731.52947557266        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  1.0279E+00  1.6399E+00  5.8162E-01  6.2498E-01  1.0662E+00  1.0053E+00  7.5164E-01  5.3553E-01  1.3173E+00  7.5423E-01
             1.0212E+00
 PARAMETER:  1.2748E-01  5.9463E-01 -4.4195E-01 -3.7003E-01  1.6406E-01  1.0528E-01 -1.8550E-01 -5.2449E-01  3.7556E-01 -1.8205E-01
             1.2102E-01
 GRADIENT:   2.4864E+00 -2.1745E+00 -1.4716E+00  3.0687E+00 -1.4193E+00 -4.5111E-01  8.7734E-01  7.9595E-01 -2.4527E-01  7.4606E-01
             2.9748E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1732.35062472493        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      780
 NPARAMETR:  1.0269E+00  1.8736E+00  4.4892E-01  4.6975E-01  1.1466E+00  1.0070E+00  6.8611E-01  2.8136E-01  1.5988E+00  7.9081E-01
             1.0197E+00
 PARAMETER:  1.2651E-01  7.2787E-01 -7.0091E-01 -6.5556E-01  2.3679E-01  1.0694E-01 -2.7672E-01 -1.1681E+00  5.6925E-01 -1.3469E-01
             1.1948E-01
 GRADIENT:   1.6477E-01  3.5255E+00 -7.5944E-01  1.7539E+00 -9.9112E-01  3.7778E-02  1.9062E-01  2.3420E-01  1.9162E-01  2.8589E-01
             5.4383E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1732.45055765859        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      958
 NPARAMETR:  1.0267E+00  1.8748E+00  4.5408E-01  4.6605E-01  1.1536E+00  1.0068E+00  6.8507E-01  1.2415E-01  1.6145E+00  7.9968E-01
             1.0191E+00
 PARAMETER:  1.2635E-01  7.2848E-01 -6.8949E-01 -6.6346E-01  2.4285E-01  1.0682E-01 -2.7824E-01 -1.9862E+00  5.7904E-01 -1.2355E-01
             1.1888E-01
 GRADIENT:  -4.5044E-02 -1.9148E+00  1.1402E-01 -9.6643E-01 -4.3127E-01  3.0516E-02  2.0461E-01  4.5118E-02  2.6341E-01  1.9187E-01
             2.0031E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1732.47651285310        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1133
 NPARAMETR:  1.0267E+00  1.8640E+00  4.5852E-01  4.7433E-01  1.1480E+00  1.0068E+00  6.8746E-01  1.2882E-02  1.5947E+00  7.9539E-01
             1.0183E+00
 PARAMETER:  1.2637E-01  7.2272E-01 -6.7976E-01 -6.4585E-01  2.3806E-01  1.0681E-01 -2.7475E-01 -4.2519E+00  5.6667E-01 -1.2892E-01
             1.1812E-01
 GRADIENT:  -5.3118E-02  3.6160E-01  7.0628E-02  9.2642E-02 -4.2419E-02  1.6522E-02 -5.4685E-02  4.8639E-04 -1.3776E-02 -8.1981E-02
            -7.2291E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1732.47662130380        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1225
 NPARAMETR:  1.0267E+00  1.8652E+00  4.5776E-01  4.7338E-01  1.1486E+00  1.0068E+00  6.8727E-01  1.0000E-02  1.5966E+00  7.9615E-01
             1.0184E+00
 PARAMETER:  1.2640E-01  7.2339E-01 -6.8140E-01 -6.4785E-01  2.3853E-01  1.0677E-01 -2.7502E-01 -4.5679E+00  5.6788E-01 -1.2797E-01
             1.1824E-01
 GRADIENT:   3.3415E-03  5.6738E-02  2.2030E-02  5.3310E-04 -5.9882E-02  1.7100E-03  3.2258E-03  0.0000E+00  1.0408E-02  7.2101E-03
            -3.6661E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1225
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.2889E-04 -3.1995E-02 -2.3156E-04  2.9455E-02 -3.7857E-02
 SE:             2.9870E-02  2.4129E-02  9.1003E-05  2.2982E-02  2.0989E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9121E-01  1.8484E-01  1.0941E-02  1.9995E-01  7.1282E-02

 ETASHRINKSD(%)  1.0000E-10  1.9165E+01  9.9695E+01  2.3009E+01  2.9684E+01
 ETASHRINKVR(%)  1.0000E-10  3.4657E+01  9.9999E+01  4.0724E+01  5.0557E+01
 EBVSHRINKSD(%)  4.2926E-01  1.8997E+01  9.9719E+01  2.4230E+01  2.9053E+01
 EBVSHRINKVR(%)  8.5667E-01  3.4385E+01  9.9999E+01  4.2590E+01  4.9665E+01
 RELATIVEINF(%)  9.9115E+01  5.2732E+00  8.9840E-05  4.7380E+00  1.1102E+01
 EPSSHRINKSD(%)  4.3751E+01
 EPSSHRINKVR(%)  6.8360E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1732.4766213038042     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -997.32579474006604     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.00
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1732.477       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.87E+00  4.58E-01  4.73E-01  1.15E+00  1.01E+00  6.87E-01  1.00E-02  1.60E+00  7.96E-01  1.02E+00
 


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
+        1.03E+03
 
 TH 2
+       -6.72E+00  4.64E+02
 
 TH 3
+        5.26E+00  2.00E+02  4.51E+02
 
 TH 4
+       -1.62E+01  3.37E+02 -3.69E+02  1.15E+03
 
 TH 5
+       -4.58E+00 -2.44E+02 -4.02E+02  3.27E+02  6.75E+02
 
 TH 6
+       -2.37E+00 -1.09E+00  1.76E+00 -2.11E+00 -1.54E+00  1.93E+02
 
 TH 7
+        9.40E-02  9.42E+00 -1.12E+01 -1.90E+01 -1.19E+01  2.68E-01  1.98E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.46E+00 -1.90E+01 -3.14E+01  6.37E+01 -3.69E+00 -3.76E-01  1.33E+01  0.00E+00  3.41E+01
 
 TH10
+        5.65E-02 -1.73E+01 -3.91E+01 -7.03E+00 -7.69E+01 -2.04E+00  2.55E+01  0.00E+00  7.32E+00  8.80E+01
 
 TH11
+       -5.32E+00 -2.20E+01 -2.72E+01  1.99E+00 -1.47E+01  3.89E+00  1.38E+01  0.00E+00  4.87E+00  2.38E+01  2.05E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.186
Stop Time:
Sat Sep 18 08:15:48 CDT 2021
