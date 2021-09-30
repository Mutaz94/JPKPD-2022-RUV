Wed Sep 29 22:03:45 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat38.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m38.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1890.05972332325        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5890E+02  3.4190E+01 -3.7362E+00  7.4875E+01  1.1009E+01  2.8186E+01 -2.0696E+01  1.0682E+01 -1.3900E+01 -2.2436E+00
            -3.8226E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1905.47968215589        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0518E+00  1.0869E+00  1.1529E+00  9.8175E-01  1.1111E+00  1.2166E+00  1.2442E+00  8.3600E-01  1.1527E+00  8.9014E-01
             1.9803E+00
 PARAMETER:  1.5055E-01  1.8333E-01  2.4227E-01  8.1579E-02  2.0538E-01  2.9604E-01  3.1847E-01 -7.9123E-02  2.4212E-01 -1.6372E-02
             7.8327E-01
 GRADIENT:   2.4549E+02  2.0540E+01 -9.0202E+00  2.8757E+01  4.0991E+00  8.3571E+01  1.4376E+01  4.6618E+00  2.3372E+01  7.3316E+00
             2.2120E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1932.92853029961        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0014E+00  8.6438E-01  8.8990E-01  1.0918E+00  9.0100E-01  1.0642E+00  1.3679E+00  5.5878E-02  8.4367E-01  1.0384E+00
             1.7333E+00
 PARAMETER:  1.0141E-01 -4.5741E-02 -1.6643E-02  1.8783E-01 -4.2540E-03  1.6222E-01  4.1330E-01 -2.7846E+00 -6.9991E-02  1.3765E-01
             6.5003E-01
 GRADIENT:   1.5124E+02  3.2953E+00 -3.2708E+01  7.8629E+01  1.6811E+01  3.4548E+01 -4.8954E+00  6.0153E-02 -1.7113E+01  2.6343E+01
             1.6712E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1948.14381431092        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  9.3725E-01  9.3487E-01  8.3727E-01  1.0102E+00  8.9381E-01  9.9696E-01  1.3251E+00  8.6817E-02  9.2864E-01  9.0496E-01
             1.4794E+00
 PARAMETER:  3.5197E-02  3.2655E-02 -7.7604E-02  1.1012E-01 -1.2265E-02  9.6960E-02  3.8148E-01 -2.3439E+00  2.5963E-02  1.3926E-04
             4.9167E-01
 GRADIENT:   6.2276E+01 -4.2629E+00 -9.2677E+00  9.9112E+00  3.2733E+00 -1.5972E+00  1.2936E+00  1.2993E-01 -6.0640E+00  8.7830E+00
             6.2073E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1956.65351561682        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      405
 NPARAMETR:  1.0034E+00  8.4719E-01  9.5480E-01  1.0969E+00  9.1419E-01  1.0384E+00  1.4987E+00  1.1746E-01  9.4101E-01  9.3972E-01
             1.3727E+00
 PARAMETER:  1.0335E-01 -6.5827E-02  5.3745E-02  1.9247E-01  1.0279E-02  1.3764E-01  5.0459E-01 -2.0417E+00  3.9196E-02  3.7825E-02
             4.1675E-01
 GRADIENT:   1.6947E+01  4.5487E+00 -1.1306E+00  8.0643E-01 -1.2588E+00 -3.0859E+00  1.6682E-01  1.3847E-01 -6.1225E-01  1.1541E+00
            -6.7460E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1957.48307919931        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      582
 NPARAMETR:  9.9007E-01  6.0375E-01  1.0476E+00  1.2494E+00  8.6727E-01  1.0488E+00  1.8577E+00  5.5082E-02  8.8013E-01  9.7205E-01
             1.3854E+00
 PARAMETER:  9.0017E-02 -4.0459E-01  1.4655E-01  3.2269E-01 -4.2407E-02  1.4765E-01  7.1935E-01 -2.7989E+00 -2.7688E-02  7.1653E-02
             4.2598E-01
 GRADIENT:  -3.2179E+00  4.0186E+00  1.2842E+00  7.2560E+00 -3.1611E+00  2.5147E+00  3.1433E-02  2.0557E-02 -5.0810E-01  2.5560E-01
             1.6065E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1957.57209785124        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      763
 NPARAMETR:  9.9091E-01  5.1931E-01  1.0713E+00  1.2933E+00  8.5388E-01  1.0411E+00  2.0280E+00  2.1805E-02  8.6475E-01  9.8581E-01
             1.3821E+00
 PARAMETER:  9.0867E-02 -5.5525E-01  1.6884E-01  3.5717E-01 -5.7961E-02  1.4025E-01  8.0703E-01 -3.7256E+00 -4.5312E-02  8.5706E-02
             4.2362E-01
 GRADIENT:   1.4533E+00 -4.0087E-01 -9.7707E-02 -2.6083E+00  7.4364E-01  2.6157E-01 -1.1771E-01  2.7656E-03  7.5338E-02 -7.1259E-03
            -1.7769E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1957.57390459373        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      944
 NPARAMETR:  9.9061E-01  5.1942E-01  1.0714E+00  1.2935E+00  8.5338E-01  1.0407E+00  2.0345E+00  1.0000E-02  8.6447E-01  9.8591E-01
             1.3825E+00
 PARAMETER:  9.0570E-02 -5.5505E-01  1.6899E-01  3.5732E-01 -5.8547E-02  1.3993E-01  8.1025E-01 -4.8187E+00 -4.5644E-02  8.5809E-02
             4.2386E-01
 GRADIENT:   8.4741E-01  4.3084E-02  5.6217E-01 -2.4437E+00 -3.0486E-01  1.3974E-01  1.9105E-01  0.0000E+00  8.2365E-02  8.4920E-02
             3.2945E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1957.57390459373        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      966
 NPARAMETR:  9.9061E-01  5.1942E-01  1.0714E+00  1.2935E+00  8.5338E-01  1.0407E+00  2.0345E+00  1.0000E-02  8.6447E-01  9.8591E-01
             1.3825E+00
 PARAMETER:  9.0570E-02 -5.5505E-01  1.6899E-01  3.5732E-01 -5.8547E-02  1.3993E-01  8.1025E-01 -4.8187E+00 -4.5644E-02  8.5809E-02
             4.2386E-01
 GRADIENT:   8.4741E-01  4.3084E-02  5.6217E-01 -2.4437E+00 -3.0486E-01  1.3974E-01  1.9105E-01  0.0000E+00  8.2365E-02  8.4920E-02
             3.2945E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      966
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.3314E-04  1.8730E-02 -2.8250E-04 -1.6781E-02 -1.2668E-02
 SE:             2.9782E-02  1.8271E-02  1.5917E-04  2.5077E-02  2.2609E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7500E-01  3.0532E-01  7.5931E-02  5.0338E-01  5.7528E-01

 ETASHRINKSD(%)  2.2637E-01  3.8788E+01  9.9467E+01  1.5988E+01  2.4256E+01
 ETASHRINKVR(%)  4.5223E-01  6.2531E+01  9.9997E+01  2.9420E+01  4.2628E+01
 EBVSHRINKSD(%)  6.0025E-01  4.2855E+01  9.9446E+01  1.4214E+01  2.1220E+01
 EBVSHRINKVR(%)  1.1969E+00  6.7345E+01  9.9997E+01  2.6408E+01  3.7937E+01
 RELATIVEINF(%)  9.7967E+01  3.4362E+00  4.5584E-04  9.1327E+00  8.4313E+00
 EPSSHRINKSD(%)  3.0784E+01
 EPSSHRINKVR(%)  5.2091E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1957.5739045937303     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1038.6353713890576     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.26
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.97
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1957.574       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.91E-01  5.19E-01  1.07E+00  1.29E+00  8.53E-01  1.04E+00  2.03E+00  1.00E-02  8.64E-01  9.86E-01  1.38E+00
 


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
+        1.03E+03
 
 TH 2
+       -1.77E+01  3.46E+02
 
 TH 3
+        2.65E+00  1.11E+02  3.10E+02
 
 TH 4
+       -2.34E+00  3.26E+02 -1.25E+02  6.62E+02
 
 TH 5
+        2.40E+00 -2.70E+02 -4.87E+02  1.20E+02  9.78E+02
 
 TH 6
+        1.04E+00 -3.82E+00  1.23E+00 -1.74E+00 -5.69E-01  1.79E+02
 
 TH 7
+        7.84E-01  2.46E+01  4.00E+00 -2.88E+00 -4.88E+00  2.09E-02  1.10E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.96E+00 -1.55E+01 -9.74E+00  1.08E+00  1.03E+01 -9.58E-01  8.74E+00  0.00E+00  1.61E+02
 
 TH10
+        1.55E+00  1.25E+01 -2.35E+01 -1.63E+01 -4.17E+01  6.02E-01  3.47E+00  0.00E+00  8.61E-01  7.66E+01
 
 TH11
+       -9.71E+00 -8.76E+00 -2.75E+01 -1.36E+01  1.01E+01  2.28E+00  2.36E+00  0.00E+00  3.82E+00  2.26E+01  2.25E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.301
Stop Time:
Wed Sep 29 22:04:09 CDT 2021
