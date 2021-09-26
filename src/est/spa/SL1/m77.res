Sat Sep 25 10:44:09 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat77.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m77.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1686.39543441834        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.5826E+01  6.9704E+00  2.0853E+01  5.7231E+00 -5.1601E+01 -1.2377E+01 -1.5861E+00  2.4042E+00  4.5460E+01  2.5315E+00
            -2.1104E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1696.92835868851        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0223E+00  1.0940E+00  9.8482E-01  9.2806E-01  1.0909E+00  1.0468E+00  1.0607E+00  9.6482E-01  6.9789E-01  9.9706E-01
             1.0895E+00
 PARAMETER:  1.2207E-01  1.8984E-01  8.4702E-02  2.5341E-02  1.8702E-01  1.4569E-01  1.5896E-01  6.4188E-02 -2.5969E-01  9.7059E-02
             1.8574E-01
 GRADIENT:   3.3477E+01 -1.2655E+01  1.0375E+01 -2.2704E+01  2.0662E+00  9.2074E+00 -1.2196E+01  2.8867E-01 -9.6870E-01 -2.4396E+00
             1.1308E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1699.58259218260        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0221E+00  1.0285E+00  6.7880E-01  9.6381E-01  8.7326E-01  1.0291E+00  1.2655E+00  6.3564E-01  5.7621E-01  8.2243E-01
             1.0572E+00
 PARAMETER:  1.2183E-01  1.2806E-01 -2.8743E-01  6.3139E-02 -3.5525E-02  1.2872E-01  3.3549E-01 -3.5313E-01 -4.5128E-01 -9.5492E-02
             1.5567E-01
 GRADIENT:   2.6964E+01  1.2623E+01 -2.2676E+01  3.6958E+01  2.7457E+01  5.5566E-01  4.0972E+00  3.5044E+00 -5.9140E+00  5.2953E+00
             4.6765E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1702.07768063318        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.0106E+00  8.7653E-01  5.7823E-01  1.0355E+00  7.1433E-01  1.0312E+00  1.3486E+00  2.8415E-01  6.2137E-01  6.6827E-01
             1.0389E+00
 PARAMETER:  1.1050E-01 -3.1781E-02 -4.4778E-01  1.3492E-01 -2.3641E-01  1.3072E-01  3.9910E-01 -1.1582E+00 -3.7583E-01 -3.0306E-01
             1.3813E-01
 GRADIENT:  -2.0013E+00  1.5382E+01 -6.9338E+00  4.2966E+01  8.2201E+00  1.3602E-01 -1.5421E+00  7.1839E-01  1.9273E+00  1.6164E+00
            -4.7214E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1702.59820143452        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  1.0167E+00  8.6619E-01  5.5665E-01  1.0261E+00  6.9602E-01  1.0349E+00  1.3629E+00  2.5603E-01  6.1956E-01  6.3724E-01
             1.0406E+00
 PARAMETER:  1.1658E-01 -4.3657E-02 -4.8581E-01  1.2576E-01 -2.6237E-01  1.3433E-01  4.0960E-01 -1.2624E+00 -3.7875E-01 -3.5060E-01
             1.3981E-01
 GRADIENT:  -3.2948E+01  3.1468E+00 -2.3496E+00  1.5745E+00  3.3141E-01 -4.3792E+00 -1.4806E+00  5.2672E-01  9.2980E-01  9.0490E-01
             5.8903E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1703.51795082154        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      524
 NPARAMETR:  1.0318E+00  7.0380E-01  5.8867E-01  1.1169E+00  6.5484E-01  1.0455E+00  1.6210E+00  1.3283E-01  5.8584E-01  6.6293E-01
             1.0391E+00
 PARAMETER:  1.3130E-01 -2.5126E-01 -4.2989E-01  2.1054E-01 -3.2337E-01  1.4451E-01  5.8306E-01 -1.9187E+00 -4.3471E-01 -3.1108E-01
             1.3838E-01
 GRADIENT:  -2.4636E-01  1.7411E+00  3.3554E+00 -6.3451E-01 -5.7326E+00  2.0926E-01 -7.0647E-01  7.2114E-02 -9.2489E-02  4.0756E-01
            -2.1002E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1703.53600249711        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      699
 NPARAMETR:  1.0316E+00  6.8133E-01  5.9432E-01  1.1287E+00  6.5266E-01  1.0450E+00  1.6685E+00  1.1431E-01  5.8247E-01  6.6693E-01
             1.0396E+00
 PARAMETER:  1.3116E-01 -2.8371E-01 -4.2033E-01  2.2104E-01 -3.2670E-01  1.4400E-01  6.1192E-01 -2.0688E+00 -4.4048E-01 -3.0506E-01
             1.3882E-01
 GRADIENT:   6.9186E-03 -4.2354E-01 -1.6134E-01 -5.9034E-03  2.7329E-01  4.3759E-02 -7.5832E-02  4.5142E-02  4.8196E-02 -2.7717E-02
            -6.2663E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1703.55874691181        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      878
 NPARAMETR:  1.0318E+00  6.9317E-01  5.8878E-01  1.1215E+00  6.5315E-01  1.0450E+00  1.6451E+00  2.7788E-02  5.8485E-01  6.6716E-01
             1.0396E+00
 PARAMETER:  1.3132E-01 -2.6648E-01 -4.2970E-01  2.1469E-01 -3.2595E-01  1.4405E-01  5.9781E-01 -3.4831E+00 -4.3640E-01 -3.0472E-01
             1.3885E-01
 GRADIENT:  -4.0741E-02 -1.9720E-01 -3.3590E-01 -2.9624E-02  3.8383E-01  1.0868E-02  4.0372E-03  2.5656E-03  8.1630E-02  6.6486E-02
             4.2928E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1703.55990915620        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1040
 NPARAMETR:  1.0318E+00  6.9422E-01  5.8884E-01  1.1211E+00  6.5346E-01  1.0450E+00  1.6435E+00  1.0000E-02  5.8479E-01  6.6726E-01
             1.0396E+00
 PARAMETER:  1.3134E-01 -2.6497E-01 -4.2960E-01  2.1430E-01 -3.2547E-01  1.4402E-01  5.9685E-01 -4.5320E+00 -4.3651E-01 -3.0457E-01
             1.3884E-01
 GRADIENT:   2.6822E-03  1.0523E-02  3.8916E-03  5.6078E-03 -1.5000E-02 -4.5073E-04  2.4391E-03  0.0000E+00  1.0582E-04  4.9750E-03
             7.9648E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1040
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.6030E-04  1.9122E-02 -5.5008E-04 -2.1433E-02  7.3520E-03
 SE:             2.9867E-02  2.3241E-02  2.4995E-04  2.3159E-02  2.1600E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9037E-01  4.1064E-01  2.7753E-02  3.5472E-01  7.3357E-01

 ETASHRINKSD(%)  1.0000E-10  2.2141E+01  9.9163E+01  2.2414E+01  2.7638E+01
 ETASHRINKVR(%)  1.0000E-10  3.9379E+01  9.9993E+01  3.9804E+01  4.7638E+01
 EBVSHRINKSD(%)  4.0933E-01  2.1458E+01  9.9210E+01  2.2530E+01  2.6903E+01
 EBVSHRINKVR(%)  8.1698E-01  3.8312E+01  9.9994E+01  3.9984E+01  4.6569E+01
 RELATIVEINF(%)  9.8874E+01  7.8976E+00  3.8475E-04  8.0209E+00  2.7165E+00
 EPSSHRINKSD(%)  4.2908E+01
 EPSSHRINKVR(%)  6.7405E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1703.5599091562030     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -968.40908259246487     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.31
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1703.560       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  6.94E-01  5.89E-01  1.12E+00  6.53E-01  1.05E+00  1.64E+00  1.00E-02  5.85E-01  6.67E-01  1.04E+00
 


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
+        9.48E+02
 
 TH 2
+       -1.02E+01  5.42E+02
 
 TH 3
+        2.65E+01  3.69E+02  2.15E+03
 
 TH 4
+       -1.67E+01  4.36E+02 -9.49E+02  1.58E+03
 
 TH 5
+       -8.88E+00 -6.41E+02 -2.43E+03  9.09E+02  3.12E+03
 
 TH 6
+       -1.34E+00 -2.08E+00  3.63E+00 -3.81E+00 -1.49E+00  1.81E+02
 
 TH 7
+        1.11E+00  3.92E+01 -2.63E+01 -2.10E+01  1.28E+01 -2.96E-01  3.24E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.22E-01 -2.47E+01 -6.75E+01 -1.21E-01  4.49E+01 -3.57E-01  1.75E+01  0.00E+00  2.37E+02
 
 TH10
+       -1.54E+00 -8.55E+00 -1.20E+02 -5.11E+01 -3.35E+01  5.91E-01  9.18E+00  0.00E+00  2.52E+01  1.41E+02
 
 TH11
+       -5.88E+00 -1.13E+01 -4.96E+01 -8.74E+00  1.86E+01  3.32E-01  4.00E+00  0.00E+00  2.19E+01  3.21E+01  2.01E+02
 
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
 #CPUT: Total CPU Time in Seconds,       17.080
Stop Time:
Sat Sep 25 10:44:32 CDT 2021
