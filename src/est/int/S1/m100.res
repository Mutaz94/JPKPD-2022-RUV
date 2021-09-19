Sat Sep 18 03:03:59 CDT 2021
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
$DATA ../../../../data/int/S1/dat100.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3735.16331322552        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.3723E+02  2.3169E+00 -1.2866E+01  2.0164E+01  5.3652E+01 -1.3451E+01 -1.7639E+01 -1.8136E+01  2.2576E+01 -4.6954E+00
            -1.9802E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3747.89796164091        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.5874E-01  9.6102E-01  1.0041E+00  1.0226E+00  9.4223E-01  1.0333E+00  1.1461E+00  1.1223E+00  8.1060E-01  9.8330E-01
             1.1128E+00
 PARAMETER:  5.7866E-02  6.0236E-02  1.0410E-01  1.2238E-01  4.0494E-02  1.3276E-01  2.3634E-01  2.1537E-01 -1.0998E-01  8.3156E-02
             2.0691E-01
 GRADIENT:   3.3962E+01  2.4321E+00 -2.2551E+00  1.5217E+01 -1.4106E+01  5.2474E+00  2.7381E+00 -1.9181E+00 -1.9770E+01  5.2282E+00
             4.6324E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3748.77921260901        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.5526E-01  9.7198E-01  1.0351E+00  1.0095E+00  9.5569E-01  1.0395E+00  1.2037E+00  1.1686E+00  9.0516E-01  9.5740E-01
             1.1081E+00
 PARAMETER:  5.4230E-02  7.1582E-02  1.3447E-01  1.0949E-01  5.4679E-02  1.3874E-01  2.8538E-01  2.5576E-01  3.5275E-04  5.6463E-02
             2.0262E-01
 GRADIENT:   2.7214E+01  4.2234E+00  9.6873E+00 -1.1802E+00 -1.7341E+01  7.9339E+00  1.2567E+01 -3.5661E-01  7.8739E+00  1.3971E+00
             4.5281E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3749.00197743796        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  9.5927E-01  9.7198E-01  1.0218E+00  1.0150E+00  9.5569E-01  1.0358E+00  1.2037E+00  1.1702E+00  8.7733E-01  9.5345E-01
             1.1081E+00
 PARAMETER:  5.8417E-02  7.1582E-02  1.2161E-01  1.1491E-01  5.4679E-02  1.3519E-01  2.8538E-01  2.5717E-01 -3.0874E-02  5.2333E-02
             2.0262E-01
 GRADIENT:  -1.3873E-01  1.0929E+00 -1.7952E-02  6.2892E-01 -1.3093E+01 -3.7731E-02  1.0212E+01  2.3577E-02  4.7405E-02  1.0887E-01
             4.4035E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3749.80210362726        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      507
 NPARAMETR:  9.5826E-01  9.7113E-01  1.0219E+00  1.0144E+00  9.6424E-01  1.0349E+00  1.1299E+00  1.1713E+00  8.7767E-01  9.5392E-01
             1.0828E+00
 PARAMETER:  5.7359E-02  7.0700E-02  1.2169E-01  1.1428E-01  6.3583E-02  1.3427E-01  2.2217E-01  2.5811E-01 -3.0482E-02  5.2821E-02
             1.7953E-01
 GRADIENT:   3.5441E+01 -2.4290E-01 -1.3915E+00  1.2339E+01  6.0832E+00  6.1550E+00  3.7800E-01 -1.4653E+00 -1.4116E+00 -2.3422E+00
            -3.9921E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3749.80883566167        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:      654             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5651E-01  9.7165E-01  1.0222E+00  1.0137E+00  9.6429E-01  1.0388E+00  1.1291E+00  1.1739E+00  8.7833E-01  9.5567E-01
             1.0830E+00
 PARAMETER:  5.5533E-02  7.1245E-02  1.2200E-01  1.1356E-01  6.3639E-02  1.3807E-01  2.2146E-01  2.6033E-01 -2.9731E-02  5.4653E-02
             1.7974E-01
 GRADIENT:   3.1715E+01 -2.0883E-01 -1.2827E+00  1.0965E+01  5.1309E+00  7.8539E+00  3.3485E-01 -1.2490E+00 -1.2173E+00 -2.0548E+00
            -3.3733E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3749.81715299956        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      811            RESET HESSIAN, TYPE II
 NPARAMETR:  9.5650E-01  9.7165E-01  1.0222E+00  1.0136E+00  9.6429E-01  1.0361E+00  1.1291E+00  1.1793E+00  8.7832E-01  9.5567E-01
             1.0835E+00
 PARAMETER:  5.5521E-02  7.1245E-02  1.2200E-01  1.1355E-01  6.3639E-02  1.3551E-01  2.2146E-01  2.6489E-01 -2.9744E-02  5.4653E-02
             1.8019E-01
 GRADIENT:   3.1611E+01 -2.2235E-01 -1.6283E+00  1.0827E+01  4.9701E+00  6.7169E+00  3.3414E-01 -8.3418E-01 -1.1454E+00 -2.0268E+00
            -2.0945E+00

0ITERATION NO.:   32    OBJECTIVE VALUE:  -3749.81715299956        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:      877
 NPARAMETR:  9.5650E-01  9.7165E-01  1.0222E+00  1.0136E+00  9.6429E-01  1.0361E+00  1.1291E+00  1.1793E+00  8.7832E-01  9.5567E-01
             1.0835E+00
 PARAMETER:  5.5521E-02  7.1245E-02  1.2200E-01  1.1355E-01  6.3639E-02  1.3551E-01  2.2146E-01  2.6489E-01 -2.9744E-02  5.4653E-02
             1.8019E-01
 GRADIENT:  -4.8792E+00 -5.8286E+05  7.8497E+01 -1.0266E+06 -5.3647E+00 -1.9178E-04 -2.4020E+00 -3.7820E+02 -4.3298E+01 -5.8286E+05
             6.4690E+05

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      877
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1723E-03 -1.9758E-02 -2.6649E-02  1.3838E-02 -3.2371E-02
 SE:             2.9878E-02  2.3819E-02  1.9644E-02  2.7218E-02  2.3497E-02
 N:                     100         100         100         100         100

 P VAL.:         9.1544E-01  4.0681E-01  1.7490E-01  6.1117E-01  1.6831E-01

 ETASHRINKSD(%)  1.0000E-10  2.0203E+01  3.4191E+01  8.8158E+00  2.1281E+01
 ETASHRINKVR(%)  1.0000E-10  3.6324E+01  5.6691E+01  1.6854E+01  3.8032E+01
 EBVSHRINKSD(%)  2.9611E-01  2.0154E+01  3.6378E+01  1.0874E+01  2.1103E+01
 EBVSHRINKVR(%)  5.9134E-01  3.6247E+01  5.9523E+01  2.0566E+01  3.7753E+01
 RELATIVEINF(%)  9.9407E+01  3.2709E+01  2.6872E+01  5.2285E+01  2.5769E+01
 EPSSHRINKSD(%)  2.1452E+01
 EPSSHRINKVR(%)  3.8303E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3749.8171529995629     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2095.7277932311522     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.96
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.26
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3749.817       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.56E-01  9.72E-01  1.02E+00  1.01E+00  9.64E-01  1.04E+00  1.13E+00  1.18E+00  8.78E-01  9.56E-01  1.08E+00
 


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
+        3.19E+09
 
 TH 2
+       -2.32E+00  3.09E+09
 
 TH 3
+       -2.00E+03 -2.40E+09  1.87E+09
 
 TH 4
+        2.65E+09 -2.61E+09  1.06E+05  2.20E+09
 
 TH 5
+        3.16E+09 -1.59E+04  2.42E+09  2.63E+09  3.13E+09
 
 TH 6
+        5.92E+00 -2.37E+03  4.44E+00 -2.00E+03 -3.19E+00  1.83E+02
 
 TH 7
+       -1.22E+09 -1.20E+09 -8.44E+00 -1.01E+09  1.21E+09  3.03E-02  4.66E+08
 
 TH 8
+       -9.75E+08  9.60E+08 -8.44E+01 -4.24E+04 -9.67E+08 -7.49E-01 -1.13E+00  2.98E+08
 
 TH 9
+        3.47E+09 -3.91E+00  9.48E+04  2.88E+09  3.44E+09  9.07E-01 -1.33E+09 -1.06E+09  3.78E+09
 
 TH10
+        3.19E+09 -2.94E+01  2.07E+05  2.65E+09 -3.16E+09 -2.41E+03 -2.76E+03 -9.76E+08  3.47E+09  6.38E+09
 
 TH11
+       -1.56E+09  9.73E+04 -1.20E+09  1.30E+09  9.81E+04  1.18E+03  5.97E+08  4.78E+08 -1.70E+09 -1.56E+09  7.64E+08
 
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
 #CPUT: Total CPU Time in Seconds,       37.352
Stop Time:
Sat Sep 18 03:04:38 CDT 2021
