Wed Sep 29 20:21:15 CDT 2021
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
$DATA ../../../../data/spa/D/dat81.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m81.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   23348.0640928312        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.0938E+02  3.9922E+02 -5.5862E+01  5.2457E+02  5.3434E+01 -2.1280E+03 -1.0774E+03 -2.5876E+01 -1.3299E+03 -2.2593E+02
            -4.4630E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -524.350510982286        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1752E+00  1.0510E+00  9.7353E-01  1.2211E+00  1.2143E+00  1.5866E+00  1.0983E+00  9.7467E-01  9.9737E-01  9.6548E-01
             1.4928E+01
 PARAMETER:  2.6141E-01  1.4973E-01  7.3171E-02  2.9972E-01  2.9417E-01  5.6158E-01  1.9378E-01  7.4347E-02  9.7363E-02  6.4866E-02
             2.8033E+00
 GRADIENT:  -3.9498E+01 -2.1160E+00 -1.1755E+00 -1.0220E+01 -6.4598E+00  2.3174E+01  1.7623E+00  2.3855E+00  8.5131E+00  2.7719E+00
             6.8569E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -533.189446742398        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.2300E+00  8.7449E-01  1.0216E+00  1.4554E+00  2.4666E+00  1.4619E+00  1.7535E+00  3.1671E-01  6.7581E-01  9.5617E-01
             1.4999E+01
 PARAMETER:  3.0705E-01 -3.4115E-02  1.2136E-01  4.7530E-01  1.0028E+00  4.7976E-01  6.6159E-01 -1.0498E+00 -2.9185E-01  5.5185E-02
             2.8080E+00
 GRADIENT:  -1.1562E+01  1.5965E+01 -2.0651E+00  2.4085E+01 -6.4022E+00 -7.9176E+00  2.3984E+00  2.4964E-01  1.0568E+00  2.1655E-02
             3.4494E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -557.493552577345        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  1.0918E+00  4.6039E-01  6.1985E-01  1.3895E+00  2.2254E+01  1.5508E+00  2.0445E+00  1.0000E-02  5.3242E-01  1.0478E+01
             1.3441E+01
 PARAMETER:  1.8781E-01 -6.7569E-01 -3.7828E-01  4.2896E-01  3.2025E+00  5.3877E-01  8.1515E-01 -6.3466E+00 -5.3032E-01  2.4493E+00
             2.6983E+00
 GRADIENT:  -4.3933E+00  1.2798E+01  2.0768E+01 -3.1060E+00  7.2658E+00  1.4425E+01  1.8676E+00  0.0000E+00  4.0282E+00 -1.4657E+00
            -1.8884E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -572.965278479198        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      311
 NPARAMETR:  1.0005E+00  1.4707E-01  4.4248E-01  1.3897E+00  1.4296E+01  1.4539E+00  1.2880E+00  1.0000E-02  1.0832E-01  9.3345E+00
             1.4011E+01
 PARAMETER:  1.0047E-01 -1.8169E+00 -7.1537E-01  4.2910E-01  2.7600E+00  4.7423E-01  3.5308E-01 -9.8532E+00 -2.1227E+00  2.3337E+00
             2.7398E+00
 GRADIENT:  -5.0361E+01  1.6234E+00  2.6753E+01 -9.3495E+00  1.1199E+01 -1.2834E+01  1.2285E-01  0.0000E+00  6.7782E-01 -1.1552E+01
             4.1530E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -576.390438298834        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:      455             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0772E+00  8.0976E-02  3.4354E-01  1.3783E+00  1.2613E+01  1.5200E+00  9.8570E-01  1.0000E-02  2.3628E-02  8.4770E+00
             1.3954E+01
 PARAMETER:  1.7433E-01 -2.4136E+00 -9.6846E-01  4.2084E-01  2.6347E+00  5.1873E-01  8.5595E-02 -1.2178E+01 -3.6453E+00  2.2374E+00
             2.7357E+00
 GRADIENT:   3.5821E+01  5.8237E-01  4.6184E+00  4.7482E+00  1.3335E+01 -7.4483E+00  2.1120E-02  0.0000E+00  3.6057E-02 -1.5152E+01
             6.7648E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -578.335876397996        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      553
 NPARAMETR:  1.0363E+00  8.0902E-02  3.4353E-01  1.3779E+00  1.2533E+01  1.5200E+00  8.9525E-01  1.0000E-02  2.0607E-02  8.5047E+00
             1.3929E+01
 PARAMETER:  1.3563E-01 -2.4145E+00 -9.6849E-01  4.2058E-01  2.6283E+00  5.1873E-01 -1.0653E-02 -1.2178E+01 -3.7821E+00  2.2406E+00
             2.7339E+00
 GRADIENT:  -1.3993E+00  1.2156E+00  2.6109E+00  1.1806E+01  1.0705E+01 -1.1921E+01  1.4563E-02  0.0000E+00  2.3939E-02 -1.7636E+01
            -9.3317E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -578.422824788269        NO. OF FUNC. EVALS.: 107
 CUMULATIVE NO. OF FUNC. EVALS.:      660
 NPARAMETR:  1.0359E+00  8.0234E-02  3.4438E-01  1.3757E+00  1.2485E+01  1.5180E+00  8.7157E-01  1.0000E-02  1.0043E-02  8.5102E+00
             1.3977E+01
 PARAMETER:  1.3523E-01 -2.4228E+00 -9.6602E-01  4.1895E-01  2.6245E+00  5.1742E-01 -3.7460E-02 -1.2178E+01 -4.5008E+00  2.2413E+00
             2.7374E+00
 GRADIENT:   4.2492E+00  1.6571E+00  5.5567E+00  1.6356E+01  8.6216E+00 -5.9545E+00  7.6956E-03  0.0000E+00  6.6878E-03 -1.2029E+01
             1.5263E+01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -578.422824788269        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  1.0358E+00  8.0176E-02  3.4446E-01  1.3755E+00  1.2482E+01  1.5178E+00  8.7070E-01  1.0000E-02  1.0026E-02  8.5082E+00
             1.3985E+01
 PARAMETER:  1.3523E-01 -2.4228E+00 -9.6602E-01  4.1895E-01  2.6245E+00  5.1742E-01 -3.7460E-02 -1.2178E+01 -4.5008E+00  2.2413E+00
             2.7374E+00
 GRADIENT:   4.9061E+02  1.4762E+01 -6.4467E+01  8.4489E+01  1.3013E+01  1.1693E+02  7.8069E-03  0.0000E+00  1.8855E-03  1.1468E+01
            -2.7578E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      725
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.8991E-03 -6.8104E-04  3.3634E-05 -5.4445E-04 -3.1252E-02
 SE:             2.8883E-02  3.3510E-04  1.0021E-04  2.5771E-04  8.4606E-03
 N:                     100         100         100         100         100

 P VAL.:         8.1121E-01  4.2117E-02  7.3715E-01  3.4630E-02  2.2090E-04

 ETASHRINKSD(%)  3.2381E+00  9.8877E+01  9.9664E+01  9.9137E+01  7.1656E+01
 ETASHRINKVR(%)  6.3714E+00  9.9987E+01  9.9999E+01  9.9993E+01  9.1966E+01
 EBVSHRINKSD(%)  6.4164E+00  9.8510E+01  9.9607E+01  9.9048E+01  7.1911E+01
 EBVSHRINKVR(%)  1.2421E+01  9.9978E+01  9.9998E+01  9.9991E+01  9.2110E+01
 RELATIVEINF(%)  3.5478E+01  1.6226E-03  9.5899E-05  2.3303E-04  2.8879E+00
 EPSSHRINKSD(%)  2.3432E+00
 EPSSHRINKVR(%)  4.6315E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -578.42282478826860     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       156.72800177546958     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.86
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.90
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -578.423       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  8.02E-02  3.44E-01  1.38E+00  1.25E+01  1.52E+00  8.72E-01  1.00E-02  1.00E-02  8.51E+00  1.40E+01
 


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
+        8.52E+04
 
 TH 2
+       -2.09E+02  4.42E+04
 
 TH 3
+        4.20E+01  4.90E+01  1.57E+04
 
 TH 4
+       -1.49E+02 -2.61E+02 -2.17E+02  5.26E+03
 
 TH 5
+        1.95E+00 -5.94E+02  1.83E+02 -1.02E+02  3.38E+00
 
 TH 6
+       -2.61E+01 -6.18E+01  8.79E+01 -6.46E+01  7.90E+01  2.76E+03
 
 TH 7
+        4.01E-01  3.95E-01 -3.64E-02  3.31E-02 -4.95E-04  5.85E-02  1.07E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.11E+00  2.11E+00 -1.18E+00  5.23E-01 -2.68E-02  5.54E-01  6.07E-01  0.00E+00  3.26E+02
 
 TH10
+       -2.26E+00  1.16E+00 -3.98E-02 -5.80E+00 -6.10E+00 -3.41E-01  3.55E-04  0.00E+00  3.19E-02  1.50E+01
 
 TH11
+       -8.85E+00 -1.67E+00  1.19E+01 -7.07E+00  2.63E+00  1.91E+00 -8.29E-04  0.00E+00  1.09E-02 -4.40E+00  5.57E+00
 
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
 #CPUT: Total CPU Time in Seconds,       19.834
Stop Time:
Wed Sep 29 20:21:36 CDT 2021
