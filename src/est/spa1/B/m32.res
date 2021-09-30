Wed Sep 29 21:02:27 CDT 2021
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
$DATA ../../../../data/spa1/B/dat32.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m32.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1867.70211838590        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4882E+02 -2.7240E+01  1.0129E+02 -6.0879E+01  4.1835E+01  8.0695E+01  5.4972E+00 -4.3144E+02 -8.9543E+01  7.5245E+00
            -3.7140E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2075.84286995057        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.5412E-01  1.1399E+00  1.0112E+00  9.2494E-01  9.8115E-01  1.0293E+00  9.0355E-01  2.8435E+00  9.2615E-01  8.0764E-01
             9.6574E-01
 PARAMETER:  5.3032E-02  2.3091E-01  1.1114E-01  2.1973E-02  8.0975E-02  1.2892E-01 -1.4228E-03  1.1451E+00  2.3283E-02 -1.1364E-01
             6.5142E-02
 GRADIENT:   2.8875E+02  6.5676E+01  6.3593E+00 -3.1064E+01 -5.9270E+01  1.0788E+02  3.7525E+00  9.9070E+00 -1.3855E+01  1.0921E+01
            -3.1505E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2085.96332836571        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  9.5163E-01  1.0474E+00  1.4933E+00  1.0167E+00  1.1037E+00  1.0112E+00  8.4034E-01  3.2126E+00  9.0846E-01  7.7392E-01
             1.0203E+00
 PARAMETER:  5.0420E-02  1.4635E-01  5.0101E-01  1.1659E-01  1.9870E-01  1.1115E-01 -7.3947E-02  1.2671E+00  3.9991E-03 -1.5628E-01
             1.2014E-01
 GRADIENT:   2.3506E+02  1.8740E+01  7.3218E+00  4.2084E+00 -9.8481E+00  8.6546E+01 -1.2664E-01  1.8448E+01 -1.4583E+01 -5.4165E+00
            -1.9815E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2092.73571877428        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      340
 NPARAMETR:  9.6529E-01  1.0555E+00  1.6767E+00  1.0234E+00  1.1611E+00  9.9370E-01  8.6019E-01  3.3987E+00  9.0573E-01  8.6909E-01
             1.0317E+00
 PARAMETER:  6.4673E-02  1.5405E-01  6.1681E-01  1.2312E-01  2.4936E-01  9.3675E-02 -5.0604E-02  1.3234E+00  9.8275E-04 -4.0306E-02
             1.3119E-01
 GRADIENT:  -1.5954E+02 -5.1281E+01  1.2405E-01 -6.0036E+01 -7.5684E+00  2.8100E+01 -7.2217E+00 -2.8546E+01 -1.5437E+01 -6.6238E+00
             1.4222E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2101.15329404185        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      496
 NPARAMETR:  1.0191E+00  1.0558E+00  1.6752E+00  1.0236E+00  1.1607E+00  8.9388E-01  9.1194E-01  3.3989E+00  9.4487E-01  9.0048E-01
             1.0315E+00
 PARAMETER:  1.1888E-01  1.5428E-01  6.1591E-01  1.2330E-01  2.4900E-01 -1.2180E-02  7.8216E-03  1.3234E+00  4.3291E-02 -4.8321E-03
             1.3100E-01
 GRADIENT:  -4.6496E+01 -4.4587E+01 -1.1022E+00 -5.6940E+01 -2.2893E+01  1.8679E+00  2.9240E-02 -2.4905E+01 -3.3490E+00 -5.9371E-01
             1.0043E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2101.63397970076        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      671
 NPARAMETR:  1.0363E+00  1.0557E+00  1.6758E+00  1.0235E+00  1.1609E+00  8.8816E-01  8.9424E-01  3.3997E+00  9.6535E-01  9.1414E-01
             1.0316E+00
 PARAMETER:  1.3569E-01  1.5422E-01  6.1627E-01  1.2325E-01  2.4915E-01 -1.8600E-02 -1.1777E-02  1.3237E+00  6.4731E-02  1.0229E-02
             1.3107E-01
 GRADIENT:   4.9790E-01 -4.7307E+01 -7.9554E-01 -5.3641E+01 -2.2380E+01  1.5541E-01  9.9059E-02 -2.5420E+01 -7.5447E-02  3.2166E-01
             1.2455E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2106.08708547558        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      849
 NPARAMETR:  1.0407E+00  1.1613E+00  1.8116E+00  9.6204E-01  1.2406E+00  8.9545E-01  9.0137E-01  3.8198E+00  9.3573E-01  1.0016E+00
             9.8960E-01
 PARAMETER:  1.3992E-01  2.4955E-01  6.9419E-01  6.1303E-02  3.1562E-01 -1.0426E-02 -3.8448E-03  1.4402E+00  3.3567E-02  1.0155E-01
             8.9541E-02
 GRADIENT:   1.0886E+01 -5.2593E+01 -1.2247E+00 -4.5396E+01  3.9728E-01  2.6064E+00 -2.4135E+00 -2.0113E+01 -7.9184E+00 -2.8231E+00
            -1.3330E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2144.76894628414        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1038
 NPARAMETR:  1.0385E+00  1.1823E+00  1.8296E+00  9.5173E-01  1.2549E+00  8.9135E-01  9.0276E-01  3.9231E+00  9.3724E-01  1.0170E+00
             9.8744E-01
 PARAMETER:  1.3775E-01  2.6746E-01  7.0409E-01  5.0528E-02  3.2703E-01 -1.5017E-02 -2.2949E-03  1.4669E+00  3.5183E-02  1.1687E-01
             8.7356E-02
 GRADIENT:   4.2247E+00 -2.9776E+01 -1.2392E+01 -2.8096E+01 -4.1287E-01  2.8975E-01  1.5883E+00  4.7806E+01 -8.1676E+00 -3.6951E+00
            -9.2232E+00

0ITERATION NO.:   36    OBJECTIVE VALUE:  -2144.76894628414        NO. OF FUNC. EVALS.:  25
 CUMULATIVE NO. OF FUNC. EVALS.:     1063
 NPARAMETR:  1.0385E+00  1.1823E+00  1.8296E+00  9.5173E-01  1.2549E+00  8.9135E-01  9.0276E-01  3.9231E+00  9.3724E-01  1.0170E+00
             9.8744E-01
 PARAMETER:  1.3775E-01  2.6746E-01  7.0409E-01  5.0528E-02  3.2703E-01 -1.5017E-02 -2.2949E-03  1.4669E+00  3.5183E-02  1.1687E-01
             8.7356E-02
 GRADIENT:   1.2317E+00 -4.8730E+01 -2.1627E+00 -4.0293E+01  1.5685E+00  1.3289E-01 -1.9698E+00 -6.4839E+01 -7.9025E+00 -2.6485E+00
            -1.5266E+01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1063
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9210E-04 -1.1010E-02 -4.3903E-02  1.9535E-02 -7.8436E-02
 SE:             2.9882E-02  1.7926E-02  1.5180E-02  2.4846E-02  1.8617E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9487E-01  5.3908E-01  3.8251E-03  4.3173E-01  2.5188E-05

 ETASHRINKSD(%)  1.0000E-10  3.9947E+01  4.9146E+01  1.6763E+01  3.7632E+01
 ETASHRINKVR(%)  1.0000E-10  6.3937E+01  7.4139E+01  3.0717E+01  6.1103E+01
 EBVSHRINKSD(%)  3.9246E-01  3.8654E+01  3.0761E+01  2.1357E+01  3.4742E+01
 EBVSHRINKVR(%)  7.8338E-01  6.2366E+01  5.2059E+01  3.8154E+01  5.7414E+01
 RELATIVEINF(%)  9.8763E+01  1.6982E+00  1.6995E+01  2.8430E+00  2.0712E+01
 EPSSHRINKSD(%)  3.5124E+01
 EPSSHRINKVR(%)  5.7911E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2144.7689462841436     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1225.8304130794709     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.98
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2144.769       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.18E+00  1.83E+00  9.52E-01  1.25E+00  8.91E-01  9.03E-01  3.92E+00  9.37E-01  1.02E+00  9.87E-01
 


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
+        4.62E+06
 
 TH 2
+       -1.11E+01  9.45E+05
 
 TH 3
+       -4.62E-01  9.68E+00  8.64E+04
 
 TH 4
+       -4.37E+00  5.24E+02 -1.55E+01  1.04E+07
 
 TH 5
+        1.12E+00 -7.65E+01 -2.40E+01  2.35E+01  5.61E+05
 
 TH 6
+       -9.51E-01 -1.93E+00  1.34E-02 -1.20E+00  6.93E-03  5.94E+06
 
 TH 7
+        3.89E-01  2.07E+01  8.09E-01  1.16E+01 -1.27E+06  5.87E+06  1.16E+07
 
 TH 8
+       -5.81E+04  5.36E+04 -1.32E+04  1.78E+05  3.97E+04 -7.30E-02 -4.20E-01  6.90E+03
 
 TH 9
+        9.63E-01  4.69E+00 -3.01E-02  4.91E+01  7.51E+00 -6.13E-01  5.21E+01  1.27E+00  1.07E+07
 
 TH10
+        3.02E+00  5.40E+00 -1.51E+00 -5.94E+00 -8.90E+01  1.19E+00  1.51E+00  2.58E+00 -1.17E+00  6.68E+06
 
 TH11
+       -8.60E+00 -2.01E+01  4.22E-01  2.30E+00 -1.17E+06  5.36E+06  1.06E+07  1.16E+00  1.29E+01  9.82E+00  9.68E+06
 
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
 #CPUT: Total CPU Time in Seconds,       26.637
Stop Time:
Wed Sep 29 21:02:55 CDT 2021
