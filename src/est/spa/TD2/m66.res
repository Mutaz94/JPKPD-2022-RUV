Sat Sep 18 14:49:14 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat66.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m66.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1640.86995370560        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4094E+02 -7.7226E+01 -4.0904E+00 -9.5693E+01 -4.6158E+00  1.2279E+01 -8.1023E+00  8.3173E+00  2.0133E+01 -1.5283E+00
            -7.7558E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1650.95542754598        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.6336E-01  1.1651E+00  1.0417E+00  9.5957E-01  1.1147E+00  9.3349E-01  1.1352E+00  9.0151E-01  7.5866E-01  1.0677E+00
             1.0509E+00
 PARAMETER:  6.2675E-02  2.5277E-01  1.4085E-01  5.8733E-02  2.0858E-01  3.1179E-02  2.2685E-01 -3.6816E-03 -1.7621E-01  1.6552E-01
             1.4960E-01
 GRADIENT:   5.5576E+01  2.2339E+00  1.0643E+01 -1.7981E+01  4.6238E+00 -9.9799E+00 -4.1332E+00 -2.7838E-01 -7.8355E+00 -7.2111E+00
             5.2412E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1652.63274543954        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.6928E-01  8.8818E-01  1.0087E+00  1.1539E+00  9.5404E-01  9.1941E-01  1.5196E+00  5.3802E-01  6.6352E-01  1.0208E+00
             1.0510E+00
 PARAMETER:  6.8795E-02 -1.8578E-02  1.0870E-01  2.4314E-01  5.2947E-02  1.5981E-02  5.1845E-01 -5.1986E-01 -3.1019E-01  1.2063E-01
             1.4979E-01
 GRADIENT:   7.4189E+01  2.7024E+01  1.0760E+01  4.6191E+01 -5.9152E+00 -1.6604E+01  3.7027E+00 -4.1934E-01 -6.9344E+00 -1.5094E+00
             7.0009E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1654.45240276807        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.4597E-01  9.4104E-01  8.1515E-01  1.0894E+00  8.7453E-01  9.4748E-01  1.3666E+00  3.1233E-01  7.4316E-01  9.2328E-01
             1.0230E+00
 PARAMETER:  4.4456E-02  3.9233E-02 -1.0439E-01  1.8565E-01 -3.4066E-02  4.6046E-02  4.1229E-01 -1.0637E+00 -1.9684E-01  2.0181E-02
             1.2275E-01
 GRADIENT:   9.0154E+00  2.4821E+00 -7.2632E+00  1.3411E+01  7.5563E+00 -4.1484E+00  1.2112E+00  6.8415E-01  1.5250E+00  1.6428E+00
             1.3426E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1654.46769381217        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  9.4284E-01  9.3585E-01  7.9638E-01  1.0860E+00  8.5858E-01  9.5551E-01  1.3652E+00  2.6040E-01  7.4085E-01  9.0255E-01
             1.0214E+00
 PARAMETER:  4.1140E-02  3.3700E-02 -1.2768E-01  1.8247E-01 -5.2476E-02  5.4487E-02  4.1127E-01 -1.2456E+00 -1.9996E-01 -2.5362E-03
             1.2120E-01
 GRADIENT:   1.2101E+00 -2.9115E-01 -3.2101E+00  2.9635E+00  2.7736E+00 -9.3361E-01  4.7738E-01  4.6130E-01  8.8708E-01  8.8364E-01
             4.9564E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1654.69964776158        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  9.4087E-01  8.8061E-01  7.6723E-01  1.1100E+00  8.1423E-01  9.6157E-01  1.4264E+00  3.6278E-02  7.2204E-01  8.6608E-01
             1.0203E+00
 PARAMETER:  3.9047E-02 -2.7145E-02 -1.6497E-01  2.0433E-01 -1.0551E-01  6.0813E-02  4.5512E-01 -3.2165E+00 -2.2567E-01 -4.3783E-02
             1.2011E-01
 GRADIENT:  -3.5302E+00 -2.8683E+00  1.0519E+00 -5.3828E+00 -2.4077E+00  1.4438E+00 -3.4336E-01  9.6619E-03 -1.2391E-01  2.7100E-01
            -2.2907E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1655.61139467917        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      532
 NPARAMETR:  9.5601E-01  6.9658E-01  8.2039E-01  1.2406E+00  7.6913E-01  9.6080E-01  1.7682E+00  1.0000E-02  6.7284E-01  8.8255E-01
             1.0197E+00
 PARAMETER:  5.5010E-02 -2.6158E-01 -9.7980E-02  3.1560E-01 -1.6250E-01  6.0011E-02  6.6996E-01 -6.8613E+00 -2.9625E-01 -2.4944E-02
             1.1954E-01
 GRADIENT:   4.1418E+00  1.1685E+01  3.2661E+00  2.4681E+01 -8.9560E+00 -1.4345E+00 -1.4451E-01  0.0000E+00 -1.3407E+00  3.4044E-01
            -1.3595E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1655.86582377916        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      708
 NPARAMETR:  9.5329E-01  6.2765E-01  8.3699E-01  1.2676E+00  7.6225E-01  9.6379E-01  1.9040E+00  1.0000E-02  6.6361E-01  8.8614E-01
             1.0220E+00
 PARAMETER:  5.2163E-02 -3.6577E-01 -7.7940E-02  3.3712E-01 -1.7147E-01  6.3116E-02  7.4398E-01 -9.0751E+00 -3.1007E-01 -2.0880E-02
             1.2180E-01
 GRADIENT:  -3.2262E-01  6.1349E-01  4.1487E-01  1.1673E+00 -6.2429E-01  1.8180E-01 -4.0827E-02  0.0000E+00 -7.2150E-02 -1.0997E-01
            -2.2879E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1655.86732892164        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      835
 NPARAMETR:  9.5331E-01  6.1908E-01  8.3845E-01  1.2718E+00  7.6073E-01  9.6322E-01  1.9239E+00  1.0000E-02  6.6212E-01  8.8709E-01
             1.0225E+00
 PARAMETER:  5.2185E-02 -3.7953E-01 -7.6198E-02  3.4041E-01 -1.7347E-01  6.2527E-02  7.5433E-01 -9.3642E+00 -3.1231E-01 -1.9808E-02
             1.2223E-01
 GRADIENT:   2.2216E-03  2.9480E-04  4.8266E-03 -2.8431E-03 -5.8636E-03 -1.3102E-03  2.1635E-04  0.0000E+00 -9.4377E-04  1.1872E-03
            -1.4140E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      835
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.2190E-04  2.2209E-02 -4.8456E-04 -2.5404E-02 -1.0162E-03
 SE:             2.9825E-02  2.1490E-02  2.1329E-04  2.2858E-02  2.3383E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8871E-01  3.0138E-01  2.3100E-02  2.6640E-01  9.6534E-01

 ETASHRINKSD(%)  8.0909E-02  2.8005E+01  9.9285E+01  2.3423E+01  2.1664E+01
 ETASHRINKVR(%)  1.6175E-01  4.8168E+01  9.9995E+01  4.1359E+01  3.8634E+01
 EBVSHRINKSD(%)  4.8319E-01  2.8874E+01  9.9324E+01  2.2270E+01  1.9298E+01
 EBVSHRINKVR(%)  9.6404E-01  4.9410E+01  9.9995E+01  3.9581E+01  3.4871E+01
 RELATIVEINF(%)  9.8421E+01  6.6753E+00  4.4027E-04  8.6231E+00  5.5714E+00
 EPSSHRINKSD(%)  4.2643E+01
 EPSSHRINKVR(%)  6.7102E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1655.8673289216392     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -920.71650235790105     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.41
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.90
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1655.867       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.53E-01  6.19E-01  8.38E-01  1.27E+00  7.61E-01  9.63E-01  1.92E+00  1.00E-02  6.62E-01  8.87E-01  1.02E+00
 


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
+        1.31E+03
 
 TH 2
+       -1.39E+01  4.19E+02
 
 TH 3
+        1.67E+01  1.76E+02  7.90E+02
 
 TH 4
+       -7.98E+00  3.80E+02 -3.21E+02  9.73E+02
 
 TH 5
+       -4.79E+00 -3.42E+02 -1.01E+03  3.46E+02  1.59E+03
 
 TH 6
+        1.10E+00 -3.95E+00  7.98E+00 -2.59E+00 -2.32E+00  2.14E+02
 
 TH 7
+        1.84E+00  3.56E+01 -2.06E+00 -1.36E+01  1.41E-01  1.40E-01  1.98E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.51E+00 -2.02E+01 -4.25E+01 -1.92E+01  4.25E+01 -1.12E+00  1.63E+01  0.00E+00  1.79E+02
 
 TH10
+        3.71E+00 -6.57E+00 -6.01E+01 -2.87E+01 -4.80E+01 -3.37E+00  3.08E+00  0.00E+00  1.31E+01  1.06E+02
 
 TH11
+       -1.01E+01 -9.40E+00 -4.46E+01 -5.81E+00  1.63E+01  1.98E+00  2.16E+00  0.00E+00  1.60E+01  3.11E+01  2.11E+02
 
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
 #CPUT: Total CPU Time in Seconds,       14.379
Stop Time:
Sat Sep 18 14:49:30 CDT 2021
