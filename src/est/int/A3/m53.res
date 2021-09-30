Wed Sep 29 00:08:44 CDT 2021
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
$DATA ../../../../data/int/A3/dat53.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m53.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -522.203193525453        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5513E+02  7.8829E+01  1.5222E+02 -3.0601E+00  1.7452E+02  3.5479E+01 -6.8858E+01 -5.0797E+01  6.5588E+00 -7.6461E+01
            -6.3934E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2750.97389217749        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.7853E-01  1.0482E+00  8.9094E-01  1.1015E+00  9.2449E-01  9.5573E-01  9.0787E-01  9.2590E-01  8.6987E-01  1.0170E+00
             2.5957E+00
 PARAMETER:  7.8293E-02  1.4704E-01 -1.5479E-02  1.9671E-01  2.1482E-02  5.4719E-02  3.3410E-03  2.3008E-02 -3.9412E-02  1.1690E-01
             1.0539E+00
 GRADIENT:   6.4071E+01  7.2254E+01 -3.1499E+00  9.1437E+01 -6.5785E+00 -5.9047E+00  5.5915E+00  1.2610E+01 -6.2373E+00 -3.2318E+00
            -2.0489E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2758.22313765919        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  9.8248E-01  1.0044E+00  8.0302E-01  1.1092E+00  8.8079E-01  9.6434E-01  8.0562E-01  2.4120E-01  8.9610E-01  1.0063E+00
             2.6003E+00
 PARAMETER:  8.2328E-02  1.0441E-01 -1.1938E-01  2.0368E-01 -2.6939E-02  6.3686E-02 -1.1614E-01 -1.3221E+00 -9.7024E-03  1.0623E-01
             1.0556E+00
 GRADIENT:   7.3339E+01  5.3643E+01 -9.5248E+00  7.5484E+01  1.6430E+01 -2.5575E+00 -2.4473E+00  1.0199E+00 -4.6989E+00 -1.3558E+01
            -1.1715E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2761.04021449708        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  9.4917E-01  7.5619E-01  6.7263E-01  1.1935E+00  6.6484E-01  9.7079E-01  9.0075E-01  1.4548E-01  8.6460E-01  9.4227E-01
             2.5870E+00
 PARAMETER:  4.7829E-02 -1.7947E-01 -2.9656E-01  2.7686E-01 -3.0821E-01  7.0351E-02 -4.5291E-03 -1.8277E+00 -4.5486E-02  4.0533E-02
             1.0505E+00
 GRADIENT:  -4.0842E+00  2.9591E+01  1.4349E+01  3.7006E+01  9.4735E+00 -7.7832E-01 -9.8851E-01  7.4042E-01 -2.3786E+00 -1.5892E-02
            -5.4106E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2762.35837462995        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      414
 NPARAMETR:  9.7328E-01  7.1059E-01  6.4625E-01  1.2124E+00  6.3731E-01  9.9323E-01  9.3130E-01  8.4143E-02  8.6717E-01  9.2861E-01
             2.6230E+00
 PARAMETER:  7.2917E-02 -2.4166E-01 -3.3656E-01  2.9258E-01 -3.5050E-01  9.3204E-02  2.8826E-02 -2.3752E+00 -4.2524E-02  2.5931E-02
             1.0643E+00
 GRADIENT:  -5.5936E-01 -3.8933E+00  3.6171E+00 -7.5228E-01 -4.4292E-01  2.8950E+00 -8.1742E-01  2.4666E-01 -3.2662E-01  2.6380E+00
             9.4628E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2767.94059623389        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      592
 NPARAMETR:  9.8526E-01  3.8154E-01  2.9760E-01  1.2716E+00  3.1633E-01  9.6475E-01  1.3890E+00  1.0000E-02  9.3761E-01  5.2479E-01
             2.5430E+00
 PARAMETER:  8.5151E-02 -8.6353E-01 -1.1120E+00  3.4027E-01 -1.0510E+00  6.4114E-02  4.2859E-01 -7.8594E+00  3.5579E-02 -5.4476E-01
             1.0334E+00
 GRADIENT:   2.4833E+01  4.7398E+01 -3.4140E+00  1.0473E+02 -3.5881E+01 -1.1375E+01 -4.5288E+00  0.0000E+00 -1.3286E+01 -2.0013E+01
            -6.1918E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2779.40923942570        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      768
 NPARAMETR:  9.8169E-01  2.6491E-01  1.9819E-01  1.1076E+00  2.3091E-01  9.9008E-01  1.4054E+00  1.0000E-02  1.1149E+00  6.3694E-01
             2.4811E+00
 PARAMETER:  8.1518E-02 -1.2284E+00 -1.5185E+00  2.0217E-01 -1.3657E+00  9.0030E-02  4.4032E-01 -1.0387E+01  2.0877E-01 -3.5107E-01
             1.0087E+00
 GRADIENT:   1.5386E+01 -2.2908E+00  7.7390E+00 -6.9128E+00 -9.0006E-01 -7.6283E-01 -5.2194E-01  0.0000E+00  4.0256E+00  5.0116E+00
             7.0856E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2779.72197261310        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      946
 NPARAMETR:  9.7542E-01  2.5833E-01  1.8939E-01  1.0996E+00  2.2402E-01  9.9230E-01  1.4204E+00  1.0000E-02  1.1180E+00  6.2171E-01
             2.4700E+00
 PARAMETER:  7.5112E-02 -1.2535E+00 -1.5639E+00  1.9492E-01 -1.3960E+00  9.2272E-02  4.5094E-01 -1.0638E+01  2.1158E-01 -3.7528E-01
             1.0042E+00
 GRADIENT:   1.2571E+00  1.2552E+00  1.1624E+00 -2.1019E-01 -3.0938E+00  2.7216E-02  2.8343E-01  0.0000E+00  1.0190E-01  1.3081E-01
             1.3397E-01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -2779.72297094326        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     1044
 NPARAMETR:  9.7444E-01  2.5810E-01  1.8924E-01  1.0997E+00  2.2393E-01  9.9222E-01  1.4170E+00  1.0000E-02  1.1174E+00  6.2099E-01
             2.4696E+00
 PARAMETER:  7.4112E-02 -1.2544E+00 -1.5648E+00  1.9506E-01 -1.3964E+00  9.2185E-02  4.4857E-01 -1.0638E+01  2.1102E-01 -3.7644E-01
             1.0041E+00
 GRADIENT:  -9.3244E-01  7.9671E-01  6.0715E-01  3.3786E-01 -2.3198E+00 -1.6342E-02 -2.4984E-01  0.0000E+00 -1.6598E-01 -1.6221E-01
            -2.4284E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1044
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9216E-03  7.3879E-03  7.3313E-05 -7.5838E-03  4.8449E-03
 SE:             2.9449E-02  2.4753E-02  2.5444E-04  2.7532E-02  2.4764E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4797E-01  7.6535E-01  7.7325E-01  7.8297E-01  8.4489E-01

 ETASHRINKSD(%)  1.3418E+00  1.7074E+01  9.9148E+01  7.7656E+00  1.7038E+01
 ETASHRINKVR(%)  2.6655E+00  3.1233E+01  9.9993E+01  1.4928E+01  3.1173E+01
 EBVSHRINKSD(%)  1.5106E+00  1.5864E+01  9.9077E+01  6.6247E+00  1.8059E+01
 EBVSHRINKVR(%)  2.9985E+00  2.9212E+01  9.9991E+01  1.2810E+01  3.2857E+01
 RELATIVEINF(%)  9.6994E+01  1.4675E+01  4.7396E-04  3.8335E+01  3.4061E+00
 EPSSHRINKSD(%)  1.8853E+01
 EPSSHRINKVR(%)  3.4152E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2779.7229709432631     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1125.6336111748524     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.27
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.68
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2779.723       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.74E-01  2.58E-01  1.89E-01  1.10E+00  2.24E-01  9.92E-01  1.42E+00  1.00E-02  1.12E+00  6.21E-01  2.47E+00
 


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
+        1.15E+03
 
 TH 2
+        2.29E+00  5.93E+03
 
 TH 3
+        3.73E+01 -1.13E+03  2.41E+04
 
 TH 4
+       -5.34E+00 -7.91E+01 -5.58E+02  6.37E+02
 
 TH 5
+       -6.51E+00 -5.02E+03 -2.59E+04 -4.62E+02  3.86E+04
 
 TH 6
+        4.13E+00 -6.55E+00  3.02E+01 -5.77E+00 -2.22E+01  1.90E+02
 
 TH 7
+       -7.53E-01  4.58E+01  3.05E+01 -1.15E+00 -6.02E+01 -1.83E-01  4.92E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.32E+00 -3.63E+00  2.15E+02 -4.60E+00  3.81E+01 -1.17E-01  3.69E-01  0.00E+00  1.15E+02
 
 TH10
+       -3.94E+00  8.21E-01 -1.90E+01  1.44E+01  8.82E+01  6.23E-01  1.62E+01  0.00E+00  3.26E+00  2.43E+02
 
 TH11
+       -1.49E+01 -1.84E+01 -1.32E+02 -9.64E+00  1.13E+02  1.85E+00  5.83E+00  0.00E+00  6.46E+00  1.57E+01  1.79E+02
 
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
 #CPUT: Total CPU Time in Seconds,       43.065
Stop Time:
Wed Sep 29 00:09:58 CDT 2021
