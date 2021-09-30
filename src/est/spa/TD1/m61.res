Wed Sep 29 18:18:43 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat61.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m61.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1714.02890731141        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8179E+02 -1.4189E+01 -1.7177E+01  2.0842E+01  2.0220E+01  5.9993E+01  1.1227E+01  1.0698E+01  2.8626E+01  1.5334E+01
             6.0077E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1718.12671495276        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0256E+00  1.0492E+00  1.0396E+00  1.0264E+00  1.0182E+00  9.7025E-01  9.2600E-01  9.1114E-01  8.3563E-01  9.0129E-01
             1.0275E+00
 PARAMETER:  1.2523E-01  1.4804E-01  1.3883E-01  1.2609E-01  1.1800E-01  6.9801E-02  2.3114E-02  6.9362E-03 -7.9573E-02 -3.9252E-03
             1.2710E-01
 GRADIENT:   1.6135E+00  2.1720E+01  5.0681E+00  2.5664E+01  1.0140E+01 -2.2344E+00 -5.8496E+00  9.6511E-01 -1.3759E+01 -6.3679E+00
             4.8556E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1720.38449178608        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0276E+00  9.6565E-01  7.3767E-01  1.0543E+00  8.3029E-01  9.8577E-01  1.1631E+00  3.8918E-01  8.1687E-01  7.2208E-01
             1.0350E+00
 PARAMETER:  1.2721E-01  6.5042E-02 -2.0426E-01  1.5291E-01 -8.5976E-02  8.5668E-02  2.5111E-01 -8.4370E-01 -1.0227E-01 -2.2562E-01
             1.3437E-01
 GRADIENT:   5.2477E-01  4.5475E+00 -2.6686E+01  4.6701E+01  4.7450E+01  2.8129E+00  3.7268E+00  8.5753E-01  2.4445E+00 -6.7921E+00
             1.4288E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1721.93259761122        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      548
 NPARAMETR:  1.0255E+00  8.9859E-01  7.6173E-01  1.0811E+00  8.0606E-01  9.8227E-01  1.2011E+00  3.0598E-01  7.8914E-01  7.8509E-01
             9.8975E-01
 PARAMETER:  1.2518E-01 -6.9325E-03 -1.7216E-01  1.7798E-01 -1.1560E-01  8.2112E-02  2.8326E-01 -1.0842E+00 -1.3681E-01 -1.4195E-01
             8.9693E-02
 GRADIENT:  -1.7514E+00  3.4462E+00 -1.0836E+00  4.8166E+00  2.3633E+00  1.6976E+00  1.1838E+00  5.2290E-01  2.4665E-01 -3.3248E-01
            -1.5726E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1722.17802582810        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      723
 NPARAMETR:  1.0252E+00  7.5471E-01  7.6749E-01  1.1619E+00  7.5262E-01  9.7723E-01  1.3773E+00  1.4499E-01  7.4105E-01  7.9188E-01
             9.9069E-01
 PARAMETER:  1.2491E-01 -1.8142E-01 -1.6463E-01  2.5003E-01 -1.8420E-01  7.6970E-02  4.2013E-01 -1.8311E+00 -1.9969E-01 -1.3334E-01
             9.0646E-02
 GRADIENT:  -3.4161E-01  1.3341E+00  1.0895E+00  1.2326E+00 -1.9772E+00  3.2290E-01  7.7628E-02  8.4447E-02 -7.9651E-03 -6.9470E-02
            -2.0631E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1722.20303425475        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      882
 NPARAMETR:  1.0265E+00  7.4398E-01  7.6976E-01  1.1657E+00  7.5190E-01  9.7652E-01  1.3895E+00  6.7675E-02  7.3877E-01  7.9495E-01
             9.9109E-01
 PARAMETER:  1.2620E-01 -1.9574E-01 -1.6168E-01  2.5333E-01 -1.8516E-01  7.6237E-02  4.2894E-01 -2.5930E+00 -2.0276E-01 -1.2947E-01
             9.1054E-02
 GRADIENT:   2.9963E+00 -1.6165E+00  5.1391E-01 -3.2037E+00  1.0367E+00  1.1712E-01 -2.2988E-01  1.3803E-02 -4.4641E-02 -7.8400E-01
            -2.1337E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1722.21200539652        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1057
 NPARAMETR:  1.0249E+00  7.4585E-01  7.6830E-01  1.1660E+00  7.5126E-01  9.7629E-01  1.3884E+00  2.5884E-02  7.3898E-01  7.9974E-01
             9.9093E-01
 PARAMETER:  1.2455E-01 -1.9323E-01 -1.6358E-01  2.5357E-01 -1.8600E-01  7.6007E-02  4.2819E-01 -3.5541E+00 -2.0248E-01 -1.2346E-01
             9.0889E-02
 GRADIENT:  -1.0048E+00 -1.3565E-01  8.4497E-02 -3.3824E-01  6.7193E-02  3.3183E-03 -1.8785E-02  2.5605E-03  4.3236E-02  1.0389E-01
            -1.6837E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1722.21634652741        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1233
 NPARAMETR:  1.0254E+00  7.4695E-01  7.6531E-01  1.1651E+00  7.4978E-01  9.7635E-01  1.3879E+00  1.0000E-02  7.3889E-01  7.9687E-01
             9.9098E-01
 PARAMETER:  1.2509E-01 -1.9175E-01 -1.6748E-01  2.5281E-01 -1.8798E-01  7.6066E-02  4.2780E-01 -5.4291E+00 -2.0261E-01 -1.2707E-01
             9.0940E-02
 GRADIENT:   1.4137E-01 -7.9121E-03  1.6646E-02  3.1548E-02  6.5507E-02  4.3889E-03  1.7713E-02  0.0000E+00  2.9402E-02  2.4280E-02
             2.1293E-02

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1722.21938366063        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:     1366
 NPARAMETR:  1.0268E+00  7.4798E-01  7.6207E-01  1.1634E+00  7.4814E-01  9.7674E-01  1.3877E+00  1.0000E-02  7.3913E-01  7.9422E-01
             9.9091E-01
 PARAMETER:  1.2646E-01 -1.9038E-01 -1.7172E-01  2.5132E-01 -1.9016E-01  7.6467E-02  4.2766E-01 -7.1564E+00 -2.0229E-01 -1.3039E-01
             9.0872E-02
 GRADIENT:   3.2680E+00 -4.7458E-01  2.3596E-01 -1.6601E+00 -3.1936E-01  1.3929E-01  1.2887E-01  0.0000E+00  1.2013E-01  1.1895E-01
             6.8616E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1366
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.7532E-04  9.5531E-03 -4.7163E-04 -1.0919E-02 -5.5890E-03
 SE:             2.9848E-02  2.0499E-02  2.1831E-04  2.4635E-02  2.3237E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8997E-01  6.4120E-01  3.0743E-02  6.5761E-01  8.0992E-01

 ETASHRINKSD(%)  5.7879E-03  3.1325E+01  9.9269E+01  1.7470E+01  2.2154E+01
 ETASHRINKVR(%)  1.1575E-02  5.2838E+01  9.9995E+01  3.1887E+01  3.9400E+01
 EBVSHRINKSD(%)  4.3292E-01  3.1945E+01  9.9308E+01  1.7293E+01  2.0563E+01
 EBVSHRINKVR(%)  8.6396E-01  5.3685E+01  9.9995E+01  3.1596E+01  3.6897E+01
 RELATIVEINF(%)  9.8545E+01  3.5170E+00  4.1678E-04  6.5211E+00  4.1799E+00
 EPSSHRINKSD(%)  4.3416E+01
 EPSSHRINKVR(%)  6.7983E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1722.2193836606280     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -987.06855709688978     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.14
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.77
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1722.219       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  7.48E-01  7.62E-01  1.16E+00  7.48E-01  9.77E-01  1.39E+00  1.00E-02  7.39E-01  7.94E-01  9.91E-01
 


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
+        1.10E+03
 
 TH 2
+       -1.13E+01  4.99E+02
 
 TH 3
+        2.27E+01  2.94E+02  9.61E+02
 
 TH 4
+       -1.18E+01  4.15E+02 -3.71E+02  1.04E+03
 
 TH 5
+       -5.75E+00 -5.34E+02 -1.26E+03  3.87E+02  1.97E+03
 
 TH 6
+       -1.07E+00 -2.75E+00  4.65E+00 -2.19E+00 -1.77E+00  2.06E+02
 
 TH 7
+        1.36E+00  3.66E+01  4.56E-01 -8.17E+00 -8.55E+00  1.38E-01  2.93E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.87E+00 -2.41E+01 -4.50E+01  1.57E+01  2.08E+01 -5.42E-01  2.00E+01  0.00E+00  1.76E+02
 
 TH10
+       -1.54E+00 -1.48E+00 -8.21E+01 -3.21E+01 -4.26E+01  6.86E-02  9.44E+00  0.00E+00  1.45E+01  1.24E+02
 
 TH11
+       -7.32E+00 -1.39E+01 -4.68E+01 -6.88E+00  1.57E+01  2.07E+00  4.27E+00  0.00E+00  1.41E+01  2.95E+01  2.18E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.981
Stop Time:
Wed Sep 29 18:19:08 CDT 2021
