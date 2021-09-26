Fri Sep 24 22:21:40 CDT 2021
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
$DATA ../../../../data/int/A3/dat45.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 RAW OUTPUT FILE (FILE): m45.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   92.0544911657577        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.0482E+02  3.7012E+01  2.7995E+02 -4.3853E+01  2.0656E+02 -3.6471E+01 -1.4879E+02 -2.9134E+02 -1.0378E+02 -1.4343E+02
            -7.0922E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2319.55450417678        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.2151E-01  1.1395E+00  7.9283E-01  1.0309E+00  9.4007E-01  1.0793E+00  9.7823E-01  1.0227E+00  8.4622E-01  8.6491E-01
             5.2059E+00
 PARAMETER:  1.8259E-02  2.3061E-01 -1.3214E-01  1.3041E-01  3.8200E-02  1.7628E-01  7.7991E-02  1.2249E-01 -6.6971E-02 -4.5133E-02
             1.7498E+00
 GRADIENT:  -1.1756E+02 -3.5587E+01 -2.5287E+01 -2.6521E+01  1.2233E+01  9.4270E+00  1.4189E+01  1.0205E+01  7.1330E+00  2.2706E+01
             8.4046E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2405.60333978514        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  8.7335E-01  5.3226E-01  3.6787E-01  1.5136E+00  4.0260E-01  1.0892E+00  1.5664E+00  5.1843E-01  1.1128E+00  4.1746E-01
             4.3815E+00
 PARAMETER: -3.5421E-02 -5.3062E-01 -9.0004E-01  5.1450E-01 -8.0982E-01  1.8544E-01  5.4877E-01 -5.5695E-01  2.0692E-01 -7.7357E-01
             1.5774E+00
 GRADIENT:  -1.9218E+02  4.7438E+01 -2.2808E+01  2.6097E+02 -1.5145E+01 -1.4024E+01  3.9246E+01  1.0006E+01 -1.2243E+00  6.9015E+00
             6.7523E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2619.53751636350        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.4106E-01  6.8424E-01  5.1942E-01  1.2247E+00  5.5636E-01  1.0390E+00  1.1511E+00  1.2877E-01  9.3842E-01  6.5139E-01
             2.8513E+00
 PARAMETER:  3.9254E-02 -2.7945E-01 -5.5505E-01  3.0266E-01 -4.8634E-01  1.3823E-01  2.4070E-01 -1.9497E+00  3.6438E-02 -3.2865E-01
             1.1478E+00
 GRADIENT:  -6.9405E+00  3.5584E+01 -2.4787E+00  4.9714E+01 -3.8455E+00 -6.5958E+00 -2.3789E+00  4.9160E-01 -7.5250E+00 -1.5024E+00
             1.9004E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2624.78301719070        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  9.4111E-01  3.8228E-01  2.9405E-01  1.2547E+00  3.1571E-01  1.0979E+00  1.3261E+00  1.6258E-02  1.0918E+00  5.1619E-01
             2.7476E+00
 PARAMETER:  3.9310E-02 -8.6160E-01 -1.1240E+00  3.2693E-01 -1.0529E+00  1.9337E-01  3.8222E-01 -4.0192E+00  1.8785E-01 -5.6128E-01
             1.1107E+00
 GRADIENT:  -4.2066E+00  2.7390E+00  1.0423E+02  5.6647E+01 -9.9846E+01  1.0477E+01 -7.0065E-01 -7.4921E-03 -4.4985E+00 -1.8210E+01
            -2.0180E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2628.39981344871        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  9.4564E-01  3.0679E-01  2.0251E-01  1.1919E+00  2.4277E-01  1.0860E+00  1.2882E+00  1.0000E-02  1.1659E+00  6.1298E-01
             2.7493E+00
 PARAMETER:  4.4104E-02 -1.0816E+00 -1.4970E+00  2.7555E-01 -1.3156E+00  1.8251E-01  3.5321E-01 -6.2196E+00  2.5348E-01 -3.8942E-01
             1.1113E+00
 GRADIENT:   1.1014E-02  1.6048E+01  2.9683E+01  4.8348E+01 -5.8329E+01  3.8822E+00  4.1684E-01  0.0000E+00 -1.9461E+01 -8.7870E+00
             3.0816E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2630.63820495161        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      513
 NPARAMETR:  9.3168E-01  3.2923E-01  2.1869E-01  1.1734E+00  2.6307E-01  1.0848E+00  1.2979E+00  1.0000E-02  1.1633E+00  6.2317E-01
             2.7528E+00
 PARAMETER:  2.9234E-02 -1.0110E+00 -1.4201E+00  2.5989E-01 -1.2353E+00  1.8142E-01  3.6078E-01 -6.0692E+00  2.5127E-01 -3.7294E-01
             1.1126E+00
 GRADIENT:  -3.0270E+01  8.7270E+00 -1.2864E+01 -5.0014E-01 -2.0285E+00  2.5165E+00  1.7764E+00  0.0000E+00 -1.1027E+01 -8.9112E-01
             3.2906E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2631.47510157503        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      688
 NPARAMETR:  9.4556E-01  3.3538E-01  2.2845E-01  1.1812E+00  2.7110E-01  1.0756E+00  1.2909E+00  1.0000E-02  1.1908E+00  6.2392E-01
             2.7166E+00
 PARAMETER:  4.4019E-02 -9.9248E-01 -1.3764E+00  2.6654E-01 -1.2053E+00  1.7285E-01  3.5536E-01 -5.8683E+00  2.7458E-01 -3.7174E-01
             1.0994E+00
 GRADIENT:  -1.4305E+00  1.3729E-01  4.3328E-01 -2.6433E-01 -8.7884E-01  1.2144E-01 -4.5885E-02  0.0000E+00 -2.4952E-02 -1.4926E-02
             6.8941E-02

0ITERATION NO.:   37    OBJECTIVE VALUE:  -2631.47595085119        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      745
 NPARAMETR:  9.4622E-01  3.3588E-01  2.2893E-01  1.1819E+00  2.7157E-01  1.0751E+00  1.2912E+00  1.0000E-02  1.1901E+00  6.2340E-01
             2.7168E+00
 PARAMETER:  4.4719E-02 -9.9101E-01 -1.3743E+00  2.6709E-01 -1.2035E+00  1.7243E-01  3.5560E-01 -5.8554E+00  2.7402E-01 -3.7256E-01
             1.0995E+00
 GRADIENT:  -1.2389E-01 -2.1287E-02  1.1715E-01  2.8613E-02 -1.7383E-01 -2.0834E-02 -2.2569E-02  0.0000E+00  7.5735E-03 -4.4018E-02
            -1.4598E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      745
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2944E-03  2.9668E-03  7.1259E-05 -8.6633E-03  1.5643E-03
 SE:             2.9319E-02  2.3772E-02  3.0525E-04  2.7685E-02  2.3237E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3762E-01  9.0068E-01  8.1542E-01  7.5434E-01  9.4633E-01

 ETASHRINKSD(%)  1.7761E+00  2.0360E+01  9.8977E+01  7.2517E+00  2.2154E+01
 ETASHRINKVR(%)  3.5206E+00  3.6575E+01  9.9990E+01  1.3977E+01  3.9400E+01
 EBVSHRINKSD(%)  1.7437E+00  1.9111E+01  9.9080E+01  6.0294E+00  2.3081E+01
 EBVSHRINKVR(%)  3.4569E+00  3.4569E+01  9.9992E+01  1.1695E+01  4.0835E+01
 RELATIVEINF(%)  9.6515E+01  1.0508E+01  5.1887E-04  6.3950E+01  2.8915E+00
 EPSSHRINKSD(%)  1.7865E+01
 EPSSHRINKVR(%)  3.2538E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2631.4759508511925     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -977.38659108278171     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.81
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.02
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2631.476       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.46E-01  3.36E-01  2.29E-01  1.18E+00  2.72E-01  1.08E+00  1.29E+00  1.00E-02  1.19E+00  6.23E-01  2.72E+00
 


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
+        1.04E+03
 
 TH 2
+        9.68E-01  3.90E+03
 
 TH 3
+        2.43E+01 -1.93E+02  1.49E+04
 
 TH 4
+       -1.39E+00 -1.55E+01 -1.95E+02  4.97E+02
 
 TH 5
+        3.47E+00 -4.32E+03 -1.54E+04 -2.00E+02  2.28E+04
 
 TH 6
+        5.16E+00 -6.91E+00  4.26E+01 -6.26E+00 -2.16E+01  1.62E+02
 
 TH 7
+       -6.83E-01  3.15E+01  1.42E+01 -6.05E-01 -5.10E+01 -8.02E-01  5.06E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.37E+00 -9.63E+00  1.77E+02 -8.35E+00  8.07E+00 -5.65E-01 -1.00E-01  0.00E+00  1.05E+02
 
 TH10
+       -5.40E+00  9.94E-01  4.94E+01  8.25E+00  7.98E+01  1.20E-01  1.29E+01  0.00E+00  2.44E+00  1.85E+02
 
 TH11
+       -1.52E+01 -1.87E+01 -1.22E+02 -7.19E+00  8.50E+01  1.76E+00  8.39E+00  0.00E+00  5.31E+00  1.94E+01  1.49E+02
 
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
 #CPUT: Total CPU Time in Seconds,       28.946
Stop Time:
Fri Sep 24 22:22:11 CDT 2021
