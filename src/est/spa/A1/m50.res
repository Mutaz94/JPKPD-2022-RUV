Sat Sep 25 08:05:37 CDT 2021
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
$DATA ../../../../data/spa/A1/dat50.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1456.57772640850        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.8205E+01 -2.1816E+01  1.9594E+00 -5.0939E+01  2.8152E+01  4.0918E+00  4.3311E+00  8.0987E+00 -5.1887E+00 -6.4538E+00
            -4.0362E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1551.71066552865        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.8673E-01  9.6959E-01  1.1424E+00  1.0485E+00  9.9690E-01  9.7706E-01  8.5651E-01  8.5364E-01  9.5166E-01  8.2123E-01
             1.7590E+00
 PARAMETER:  8.6644E-02  6.9119E-02  2.3317E-01  1.4740E-01  9.6894E-02  7.6798E-02 -5.4890E-02 -5.8247E-02  5.0454E-02 -9.6946E-02
             6.6472E-01
 GRADIENT:  -1.0253E+01  4.0502E+00  1.8423E+01 -2.8522E+01 -2.3909E+01 -2.8320E+00  5.1785E+00  3.3267E+00 -6.9362E+00 -1.5622E+00
             2.2908E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1554.70610042290        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.9403E-01  8.0150E-01  9.1580E-01  1.1706E+00  8.3951E-01  9.8285E-01  4.1376E-01  4.3330E-01  1.0462E+00  8.6132E-01
             1.7095E+00
 PARAMETER:  9.4014E-02 -1.2126E-01  1.2044E-02  2.5748E-01 -7.4933E-02  8.2699E-02 -7.8246E-01 -7.3631E-01  1.4512E-01 -4.9291E-02
             6.3621E-01
 GRADIENT:   4.1244E+00  1.4345E+01 -8.9645E+00  3.8323E+01 -8.0831E+00 -1.3477E+00 -7.2594E-02  1.7191E+00  1.4108E+01  5.0009E+00
             4.7886E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1557.33849825186        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.9194E-01  6.2620E-01  1.0017E+00  1.2635E+00  8.1484E-01  9.7864E-01  4.6712E-01  1.4662E-01  9.0678E-01  8.7741E-01
             1.7130E+00
 PARAMETER:  9.1910E-02 -3.6808E-01  1.0173E-01  3.3387E-01 -1.0476E-01  7.8412E-02 -6.6116E-01 -1.8199E+00  2.1490E-03 -3.0785E-02
             6.3827E-01
 GRADIENT:   5.2586E+00  1.1842E+01  9.7701E+00  1.4046E+01 -2.2200E+01 -1.8198E+00 -1.1685E-01  1.4703E-01 -9.5002E-01  2.8837E+00
             3.6338E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1558.70922509100        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  9.8703E-01  3.9588E-01  9.6262E-01  1.3908E+00  7.3404E-01  9.8012E-01  3.0545E-01  1.3996E-02  8.2341E-01  8.4609E-01
             1.6832E+00
 PARAMETER:  8.6943E-02 -8.2665E-01  6.1907E-02  4.2991E-01 -2.0920E-01  7.9924E-02 -1.0860E+00 -4.1690E+00 -9.4301E-02 -6.7133E-02
             6.2068E-01
 GRADIENT:   1.5649E+00  4.7331E+00  1.0161E+00  1.5363E+01 -4.0658E+00 -7.1840E-01 -3.9112E-02  1.5211E-03 -2.4971E+00 -4.1166E-01
            -5.7097E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1559.48465769279        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  9.8493E-01  1.8933E-01  8.7634E-01  1.4843E+00  6.3608E-01  9.7243E-01  4.8538E-02  1.0000E-02  7.6508E-01  8.1533E-01
             1.6732E+00
 PARAMETER:  8.4810E-02 -1.5642E+00 -3.2006E-02  4.9492E-01 -3.5244E-01  7.2043E-02 -2.9254E+00 -8.4455E+00 -1.6777E-01 -1.0416E-01
             6.1472E-01
 GRADIENT:   5.2119E+00  1.1346E+00  8.8440E+00  1.0004E+00 -1.3291E+01 -2.8194E+00  4.3143E-05  0.0000E+00 -6.4285E-01  9.1783E-01
             3.3981E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1560.02641132097        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  9.8889E-01  1.4120E-01  9.2728E-01  1.5343E+00  6.5523E-01  9.8758E-01  2.5993E-02  1.0000E-02  7.5313E-01  8.3356E-01
             1.6649E+00
 PARAMETER:  8.8828E-02 -1.8576E+00  2.4495E-02  5.2808E-01 -3.2277E-01  8.7500E-02 -3.5499E+00 -1.0288E+01 -1.8352E-01 -8.2055E-02
             6.0975E-01
 GRADIENT:   3.3874E+00  1.2954E+00  2.4273E+00  9.0181E+00 -5.4255E+00  1.9666E+00 -1.1829E-05  0.0000E+00  2.7645E-01  2.8574E-01
            -6.9908E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1560.41027124666        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      694
 NPARAMETR:  9.8518E-01  4.9144E-02  9.0546E-01  1.5775E+00  6.2638E-01  9.8139E-01  1.0000E-02  1.0000E-02  7.2767E-01  8.2478E-01
             1.6601E+00
 PARAMETER:  8.5074E-02 -2.9130E+00  6.9286E-04  5.5586E-01 -3.6779E-01  8.1211E-02 -6.7737E+00 -1.7220E+01 -2.1791E-01 -9.2637E-02
             6.0690E-01
 GRADIENT:   4.9953E-01  2.5720E-01  2.3213E-01  2.8164E+00 -1.4046E+00  1.9309E-01  0.0000E+00  0.0000E+00 -4.7134E-01  2.4483E-02
            -1.2420E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1560.54344231963        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      871
 NPARAMETR:  9.8373E-01  1.0672E-02  9.0667E-01  1.5981E+00  6.1938E-01  9.8008E-01  1.0000E-02  1.0000E-02  7.1825E-01  8.2583E-01
             1.6582E+00
 PARAMETER:  8.3595E-02 -4.4402E+00  2.0208E-03  5.6883E-01 -3.7903E-01  7.9874E-02 -1.1755E+01 -2.7520E+01 -2.3093E-01 -9.1370E-02
             6.0573E-01
 GRADIENT:  -3.0393E-01  5.0866E-02  5.8632E-01  2.7856E+00 -1.3795E+00 -2.2441E-02  0.0000E+00  0.0000E+00 -3.9554E-01  1.0998E-01
            -1.6210E-01

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1560.54733918999        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      998
 NPARAMETR:  9.8381E-01  1.0000E-02  9.0782E-01  1.5976E+00  6.2018E-01  9.8007E-01  1.0000E-02  1.0000E-02  7.1889E-01  8.2571E-01
             1.6587E+00
 PARAMETER:  8.3677E-02 -4.5544E+00  3.2909E-03  5.6851E-01 -3.7774E-01  7.9872E-02 -1.2143E+01 -2.8302E+01 -2.3004E-01 -9.1513E-02
             6.0606E-01
 GRADIENT:   6.5401E-03  0.0000E+00 -4.4695E-05 -1.2812E-02  3.0799E-03 -1.2025E-03  0.0000E+00  0.0000E+00  3.1130E-03 -2.8847E-03
             1.6247E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      998
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.4882E-05 -2.8540E-06 -4.3924E-05 -8.3070E-03 -1.9512E-02
 SE:             2.9569E-02  1.7812E-06  1.9253E-04  2.8155E-02  2.3117E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9906E-01  1.0909E-01  8.1954E-01  7.6796E-01  3.9863E-01

 ETASHRINKSD(%)  9.3943E-01  9.9994E+01  9.9355E+01  5.6771E+00  2.2555E+01
 ETASHRINKVR(%)  1.8700E+00  1.0000E+02  9.9996E+01  1.1032E+01  4.0023E+01
 EBVSHRINKSD(%)  1.1081E+00  9.9994E+01  9.9318E+01  5.5574E+00  2.1748E+01
 EBVSHRINKVR(%)  2.2039E+00  1.0000E+02  9.9995E+01  1.0806E+01  3.8766E+01
 RELATIVEINF(%)  9.1471E+01  1.2625E-08  2.4728E-04  5.6618E+00  2.0346E+00
 EPSSHRINKSD(%)  3.7064E+01
 EPSSHRINKVR(%)  6.0390E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1560.5473391899898     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -825.39651262625159     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.96
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.53
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1560.547       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.84E-01  1.00E-02  9.08E-01  1.60E+00  6.20E-01  9.80E-01  1.00E-02  1.00E-02  7.19E-01  8.26E-01  1.66E+00
 


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
+        1.16E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        6.53E+00  0.00E+00  6.44E+02
 
 TH 4
+       -2.38E+01  0.00E+00 -8.50E+01  7.53E+02
 
 TH 5
+        1.42E+01  0.00E+00 -1.22E+03 -1.06E+02  2.59E+03
 
 TH 6
+       -3.55E-01  0.00E+00  1.34E+01 -8.27E+00 -9.53E+00  2.06E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.29E+01  0.00E+00  1.55E+01 -9.05E+00 -9.78E+00 -8.60E-02  0.00E+00  0.00E+00  3.08E+02
 
 TH10
+       -7.72E+00  0.00E+00  1.30E+01 -4.43E-01 -8.45E+01 -7.05E-01  0.00E+00  0.00E+00  8.77E-01  1.18E+02
 
 TH11
+       -1.12E+01  0.00E+00 -1.73E+01 -9.55E+00  6.37E+00  2.02E+00  0.00E+00  0.00E+00  1.31E+01  2.87E+01  9.25E+01
 
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
 #CPUT: Total CPU Time in Seconds,       15.549
Stop Time:
Sat Sep 25 08:05:54 CDT 2021
