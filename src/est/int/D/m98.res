Sat Sep 18 07:50:56 CDT 2021
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
$DATA ../../../../data/int/D/dat98.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m98.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   61984.5105212823        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0174E+03  9.3575E+02  1.7626E+01  9.3184E+02 -1.0016E+02 -3.6531E+03 -2.2278E+03 -3.9010E+01 -3.0468E+03 -6.6672E+02
            -1.2059E+05

0ITERATION NO.:    5    OBJECTIVE VALUE:  -265.822261933882        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0748E+00  2.2818E+00  8.6303E-01  1.5936E+00  9.5718E-01  5.9947E+00  5.0145E+00  9.8607E-01  1.7637E+00  1.1876E+00
             1.2994E+01
 PARAMETER:  1.7216E-01  9.2496E-01 -4.7304E-02  5.6601E-01  5.6233E-02  1.8909E+00  1.7123E+00  8.5972E-02  6.6743E-01  2.7192E-01
             2.6645E+00
 GRADIENT:  -7.5373E+00  5.1166E+01 -3.8082E+01  1.4689E+02 -1.9886E+01  1.5778E+02  9.9163E+00  4.0505E+00 -6.8209E+00  2.1573E+01
            -1.5409E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -336.903161917624        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  8.5278E-01  9.2908E+00  6.6908E+00  1.8095E+00  1.9445E+00  3.2558E+00  1.2199E+01  7.9661E-01  1.6529E+00  8.3958E-01
             1.3528E+01
 PARAMETER: -5.9259E-02  2.3290E+00  2.0007E+00  6.9307E-01  7.6500E-01  1.2804E+00  2.6014E+00 -1.2739E-01  6.0254E-01 -7.4853E-02
             2.7047E+00
 GRADIENT:  -5.7347E+01  4.2833E+01 -2.1841E+01  5.8756E+01  1.5766E+01  7.0114E+01 -2.3697E+00  2.4362E-01  1.2046E+01  9.8806E+00
            -2.4552E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -456.336250521702        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.3282E+00  1.7347E+00  1.5486E+01  1.4358E+00  2.3530E+00  2.1579E+00  6.9101E+00  8.7927E-01  1.3478E+00  3.7520E-01
             1.4001E+01
 PARAMETER:  3.8383E-01  6.5086E-01  2.8399E+00  4.6174E-01  9.5571E-01  8.6914E-01  2.0330E+00 -2.8666E-02  3.9848E-01 -8.8029E-01
             2.7391E+00
 GRADIENT:   6.3250E+01  2.3775E+01 -3.3634E+00  1.7629E+01  1.3884E+01 -1.2980E+01  3.2199E+01  4.3111E-02  8.1440E+00  1.8089E+00
            -2.6667E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -476.173784299633        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0296E+00  1.3422E+00  1.8605E+01  1.2083E+00  2.3608E+00  1.9410E+00  5.7780E+00  9.0185E-01  8.8271E-01  1.8327E-01
             1.4113E+01
 PARAMETER:  1.2913E-01  3.9431E-01  3.0234E+00  2.8918E-01  9.5898E-01  7.6318E-01  1.8541E+00 -3.3080E-03 -2.4756E-02 -1.5968E+00
             2.7471E+00
 GRADIENT:  -1.4908E+01  3.2979E+00 -9.6345E-01 -5.8772E-01  3.0376E+00  9.6217E+00  8.1817E+00  3.1951E-02  2.3772E+00  3.9321E-01
             2.6716E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -476.935065669114        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      399
 NPARAMETR:  1.0560E+00  1.2681E+00  2.0811E+01  1.1682E+00  2.3643E+00  1.8842E+00  5.6619E+00  8.5011E-01  7.4391E-01  1.7279E-01
             1.3892E+01
 PARAMETER:  1.5453E-01  3.3748E-01  3.1355E+00  2.5549E-01  9.6047E-01  7.3350E-01  1.8338E+00 -6.2388E-02 -1.9583E-01 -1.6557E+00
             2.7313E+00
 GRADIENT:  -5.7962E-01 -1.2945E+00 -4.9827E-01 -1.2373E+00  1.3207E+00 -5.0357E-01 -2.8636E+00  2.1342E-02  6.5983E-02  3.3218E-01
            -6.0842E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -477.255866328563        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      576
 NPARAMETR:  1.0614E+00  1.1558E+00  3.2362E+01  1.2461E+00  2.4143E+00  1.8920E+00  5.9363E+00  7.9091E-01  8.8711E-01  5.8753E-02
             1.3958E+01
 PARAMETER:  1.5961E-01  2.4480E-01  3.5770E+00  3.1998E-01  9.8142E-01  7.3763E-01  1.8811E+00 -1.3458E-01 -1.9782E-02 -2.7344E+00
             2.7360E+00
 GRADIENT:   5.6230E-01 -1.2003E+00 -3.8277E-01 -1.8258E+00  2.2318E+00  6.5004E-01 -7.5001E-01  8.0580E-03  6.9556E-01  3.8333E-02
             6.4748E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -477.447927799538        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      751
 NPARAMETR:  1.0593E+00  1.1153E+00  1.7479E+02  1.2900E+00  2.4925E+00  1.8904E+00  6.0828E+00  5.1635E-01  9.4234E-01  1.0000E-02
             1.3952E+01
 PARAMETER:  1.5765E-01  2.0911E-01  5.2636E+00  3.5461E-01  1.0133E+00  7.3680E-01  1.9055E+00 -5.6097E-01  4.0608E-02 -6.9021E+00
             2.7357E+00
 GRADIENT:  -6.1575E-01 -1.2891E-01 -3.6316E-02  5.6949E-01  1.0935E+00 -3.0732E-01 -1.2358E-01  1.1487E-04  1.4096E-03  0.0000E+00
            -1.2578E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -477.474539189938        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      927
 NPARAMETR:  1.0609E+00  1.1339E+00  2.4216E+03  1.2800E+00  2.5016E+00  1.8921E+00  6.0528E+00  2.5508E-01  9.2993E-01  1.0000E-02
             1.3964E+01
 PARAMETER:  1.5909E-01  2.2566E-01  7.8922E+00  3.4689E-01  1.0169E+00  7.3769E-01  1.9005E+00 -1.2662E+00  2.7355E-02 -1.3720E+01
             2.7365E+00
 GRADIENT:  -2.9766E-02 -1.0736E-02 -1.1962E-03 -3.6983E-02  2.8959E-02 -1.8659E-02  3.3737E-02  6.9752E-07  3.8287E-02  0.0000E+00
            -1.7389E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -477.475662394881        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1097
 NPARAMETR:  1.0612E+00  1.1384E+00  3.6427E+04  1.2778E+00  2.5019E+00  1.8923E+00  6.0456E+00  1.2283E-01  9.2492E-01  1.0000E-02
             1.3968E+01
 PARAMETER:  1.5947E-01  2.2971E-01  1.0603E+01  3.4483E-01  1.0178E+00  7.3766E-01  1.8992E+00 -1.9990E+00  2.1914E-02 -2.0749E+01
             2.7366E+00
 GRADIENT:   1.8404E-01  3.0644E-02 -5.2673E-07 -9.0624E-02  9.3265E-02 -2.2495E-02 -1.3875E-02 -1.9839E-02 -1.8318E-01  0.0000E+00
            -9.6989E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1097
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3729E-02  2.8163E-02 -5.1046E-10 -6.4557E-02 -3.0146E-07
 SE:             2.7681E-02  2.4420E-02  3.7356E-08  1.1948E-02  8.8481E-05
 N:                     100         100         100         100         100

 P VAL.:         3.9130E-01  2.4879E-01  9.8910E-01  6.5613E-08  9.9728E-01

 ETASHRINKSD(%)  7.2666E+00  1.8191E+01  1.0000E+02  5.9973E+01  9.9704E+01
 ETASHRINKVR(%)  1.4005E+01  3.3072E+01  1.0000E+02  8.3978E+01  9.9999E+01
 EBVSHRINKSD(%)  9.6640E+00  1.2390E+01  1.0000E+02  6.5344E+01  9.9626E+01
 EBVSHRINKVR(%)  1.8394E+01  2.3245E+01  1.0000E+02  8.7989E+01  9.9999E+01
 RELATIVEINF(%)  8.1303E+01  4.0012E+01  0.0000E+00  6.1898E+00  2.9992E-04
 EPSSHRINKSD(%)  3.2614E+00
 EPSSHRINKVR(%)  6.4165E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -477.47566239488089     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       1176.6136973735299     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    33.51
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.10
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -477.476       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.14E+00  3.64E+04  1.28E+00  2.50E+00  1.89E+00  6.04E+00  1.23E-01  9.25E-01  1.00E-02  1.40E+01
 


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
+        6.43E+02
 
 TH 2
+       -2.21E+01  2.58E+02
 
 TH 3
+        6.70E-05 -3.03E-05  7.20E-11
 
 TH 4
+       -3.87E+01 -1.26E+01  1.15E-05  2.36E+02
 
 TH 5
+        1.87E+01  5.64E+00  1.13E-05 -1.88E+01  2.52E+01
 
 TH 6
+       -2.15E+01  4.93E+00 -1.86E-05  7.17E+00 -1.60E+00  5.94E+01
 
 TH 7
+       -5.40E-02  3.42E+00  2.57E-06 -9.89E+00  4.95E-01 -7.03E-01  3.20E+00
 
 TH 8
+       -4.92E+01 -2.03E+02 -9.54E-05  3.11E+01 -1.39E+01 -3.83E+01 -2.99E-01  2.33E+02
 
 TH 9
+       -3.72E+02 -4.08E+02 -1.63E-04 -9.84E+01 -2.39E+01 -3.81E+01 -6.71E-01  1.42E+02  1.91E+03
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.98E+00 -1.22E+00 -7.56E-07 -8.07E+00  8.53E-01 -3.33E-02  2.95E-01 -1.82E+00  1.20E+00  0.00E+00  4.25E+00
 
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
 #CPUT: Total CPU Time in Seconds,       50.758
Stop Time:
Sat Sep 18 07:51:48 CDT 2021
