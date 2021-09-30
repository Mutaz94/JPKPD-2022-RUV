Wed Sep 29 14:18:57 CDT 2021
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
$DATA ../../../../data/spa/S1/dat49.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1586.83295678105        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2190E+02 -1.4314E+01 -1.1426E+01  2.2394E+01  2.2288E+01  3.8769E+01 -3.3992E+01  1.7050E+00 -1.8886E+01 -1.0549E+01
            -9.1549E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1598.70216534420        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.3857E-01  1.1036E+00  1.0395E+00  9.3512E-01  1.0940E+00  1.0657E+00  1.3452E+00  9.8561E-01  1.1176E+00  1.0641E+00
             1.0163E+00
 PARAMETER:  3.6605E-02  1.9859E-01  1.3872E-01  3.2920E-02  1.8986E-01  1.6360E-01  3.9653E-01  8.5507E-02  2.1117E-01  1.6215E-01
             1.1621E-01
 GRADIENT:  -2.9419E+01 -3.5731E+01 -5.6280E+00 -2.6607E+01  2.0582E+01  4.4350E+00 -1.9440E+00  1.2590E+00  1.0486E+01 -2.7595E+00
            -1.2193E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1599.69111228655        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.5268E-01  1.2402E+00  9.4859E-01  8.6643E-01  1.1075E+00  1.0742E+00  1.3472E+00  9.3758E-01  1.1201E+00  1.1285E+00
             1.0081E+00
 PARAMETER:  5.1519E-02  3.1524E-01  4.7225E-02 -4.3369E-02  2.0209E-01  1.7159E-01  3.9800E-01  3.5545E-02  2.1338E-01  2.2087E-01
             1.0804E-01
 GRADIENT:  -2.1163E+00 -1.1660E+01 -2.0844E+00 -1.3146E+01  2.6226E+00  7.7118E+00  9.8024E+00  2.1678E+00  6.2850E+00  6.3806E+00
            -2.5880E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1601.28278952694        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.5551E-01  1.2300E+00  7.2676E-01  8.7214E-01  9.6719E-01  1.0543E+00  1.3281E+00  4.2122E-01  1.0543E+00  9.5236E-01
             1.0154E+00
 PARAMETER:  5.4489E-02  3.0701E-01 -2.1917E-01 -3.6808E-02  6.6643E-02  1.5292E-01  3.8376E-01 -7.6461E-01  1.5286E-01  5.1192E-02
             1.1524E-01
 GRADIENT:   4.4956E-02 -1.0349E+00 -4.9594E+00  4.1811E+00  1.3954E+00 -5.2902E-01  1.4977E+00  9.8247E-01  1.7549E+00  1.8415E+00
             1.8479E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1601.62624870057        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  9.5588E-01  1.3150E+00  6.9141E-01  8.1968E-01  9.9414E-01  1.0563E+00  1.2602E+00  1.4764E-01  1.1004E+00  9.6424E-01
             1.0138E+00
 PARAMETER:  5.4874E-02  3.7381E-01 -2.6902E-01 -9.8844E-02  9.4122E-02  1.5481E-01  3.3126E-01 -1.8130E+00  1.9563E-01  6.3582E-02
             1.1375E-01
 GRADIENT:   2.6656E-01  2.7824E+00  7.5323E-01  3.5690E+00  7.7042E-01  2.0206E-01  3.0932E-01  7.0141E-02  5.7810E-01 -1.6235E-01
             1.9827E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1601.64014648278        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  9.5564E-01  1.3145E+00  6.8587E-01  8.1656E-01  9.9156E-01  1.0549E+00  1.2579E+00  1.4831E-01  1.0957E+00  9.6178E-01
             1.0129E+00
 PARAMETER:  5.4623E-02  3.7346E-01 -2.7707E-01 -1.0265E-01  9.1526E-02  1.5341E-01  3.2945E-01 -1.8084E+00  1.9136E-01  6.1027E-02
             1.1281E-01
 GRADIENT:  -2.8104E-01 -6.4300E-03  2.7528E-01  6.2564E-01  5.7680E-01 -3.8039E-01  9.8137E-02  8.3311E-02 -5.3335E-02  1.5520E-01
            -1.4187E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1601.65133662541        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1042
 NPARAMETR:  9.5580E-01  1.3175E+00  6.7822E-01  8.1357E-01  9.8765E-01  1.0562E+00  1.2546E+00  1.3741E-01  1.0960E+00  9.5312E-01
             1.0126E+00
 PARAMETER:  5.4793E-02  3.7576E-01 -2.8828E-01 -1.0633E-01  8.7576E-02  1.5470E-01  3.2681E-01 -1.8848E+00  1.9165E-01  5.1986E-02
             1.1257E-01
 GRADIENT:   4.0273E+02  2.1409E+02  5.0448E+00  5.3454E+01  7.5510E+00  1.0836E+02  2.6500E+01  1.3189E-01  1.1869E+01  5.8598E-01
             8.4713E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1601.65408982625        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1231             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5677E-01  1.3174E+00  6.7815E-01  8.1360E-01  9.8766E-01  1.0587E+00  1.2556E+00  1.3613E-01  1.0968E+00  9.5377E-01
             1.0127E+00
 PARAMETER:  5.5805E-02  3.7563E-01 -2.8839E-01 -1.0629E-01  8.7578E-02  1.5707E-01  3.2760E-01 -1.8941E+00  1.9236E-01  5.2666E-02
             1.1263E-01
 GRADIENT:   4.0467E+02  2.1384E+02  4.8888E+00  5.3434E+01  7.6177E+00  1.1070E+02  2.6754E+01  1.3104E-01  1.2087E+01  7.1702E-01
             9.2538E-01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1601.65408982625        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     1296
 NPARAMETR:  9.5677E-01  1.3174E+00  6.7815E-01  8.1360E-01  9.8766E-01  1.0587E+00  1.2556E+00  1.3613E-01  1.0968E+00  9.5377E-01
             1.0127E+00
 PARAMETER:  5.5805E-02  3.7563E-01 -2.8839E-01 -1.0629E-01  8.7578E-02  1.5707E-01  3.2760E-01 -1.8941E+00  1.9236E-01  5.2666E-02
             1.1263E-01
 GRADIENT:   3.0744E-01  4.1970E-01  1.2298E+05  7.1708E-02  6.7632E-03 -5.4863E-04  1.6632E-02 -2.0212E+04  1.2165E-02  2.2051E-02
             2.4791E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1296
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.5658E-04 -1.0455E-02 -5.6411E-03  7.1353E-03 -2.3539E-02
 SE:             2.9837E-02  2.4485E-02  1.9905E-03  2.2747E-02  2.1812E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9581E-01  6.6939E-01  4.5964E-03  7.5376E-01  2.8051E-01

 ETASHRINKSD(%)  4.1019E-02  1.7972E+01  9.3332E+01  2.3796E+01  2.6926E+01
 ETASHRINKVR(%)  8.2021E-02  3.2714E+01  9.9555E+01  4.1929E+01  4.6602E+01
 EBVSHRINKSD(%)  4.0178E-01  1.7367E+01  9.4263E+01  2.5023E+01  2.5574E+01
 EBVSHRINKVR(%)  8.0194E-01  3.1718E+01  9.9671E+01  4.3784E+01  4.4607E+01
 RELATIVEINF(%)  9.9068E+01  6.2306E+00  4.6463E-02  4.7660E+00  8.5283E+00
 EPSSHRINKSD(%)  4.4228E+01
 EPSSHRINKVR(%)  6.8895E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1601.6540898262483     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -866.50326326251013     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.55
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.91
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1601.654       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.57E-01  1.32E+00  6.78E-01  8.14E-01  9.88E-01  1.06E+00  1.26E+00  1.36E-01  1.10E+00  9.54E-01  1.01E+00
 


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
+        1.07E+03
 
 TH 2
+       -3.66E+00  2.88E+02
 
 TH 3
+       -1.13E+06  5.62E+04  2.50E+07
 
 TH 4
+       -9.69E+00  2.29E+02  8.58E+05  8.02E+02
 
 TH 5
+       -3.31E+00 -1.86E+02  4.05E+04  3.36E+02  6.94E+02
 
 TH 6
+       -3.86E-02 -6.45E-01 -5.84E+05 -3.12E+00 -1.35E+00  1.75E+02
 
 TH 7
+        7.76E-01  2.04E+01  2.21E+05 -1.55E+01  5.26E+00 -2.21E-01  6.14E+01
 
 TH 8
+        3.03E+05 -1.21E+04 -1.80E+06 -2.32E+05  4.77E+03  1.50E+05 -3.99E+04  1.48E+07
 
 TH 9
+        1.23E+00 -1.66E+01  2.24E+07  4.35E+01  6.63E-01 -3.96E-01  1.25E+01 -6.19E+05  1.91E+07
 
 TH10
+       -6.39E-01 -1.11E+01  1.39E+05 -1.46E+01 -6.14E+01  3.35E-04  1.05E+01 -2.45E+04  1.08E+01  7.37E+01
 
 TH11
+       -7.12E+00 -1.29E+01 -1.62E+06  2.53E+00 -1.72E+00  2.30E+00  5.37E+00  4.37E+05  6.81E+00  1.64E+01  2.10E+02
 
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
 #CPUT: Total CPU Time in Seconds,       29.513
Stop Time:
Wed Sep 29 14:19:28 CDT 2021
