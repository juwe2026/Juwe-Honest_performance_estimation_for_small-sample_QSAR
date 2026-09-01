const fs=require('fs'); const d=require('docx');
const {Document,Packer,Paragraph,TextRun,HeadingLevel,AlignmentType,ImageRun,Table,TableRow,TableCell,WidthType,ShadingType,PageBreak}=d;
const FONT="Calibri";
function T(t,o){ return new TextRun(Object.assign({text:t,font:FONT,size:22},o||{})); }
function I(t,o){ return new TextRun(Object.assign({text:t,italics:true,font:FONT,size:22},o||{})); }
function rich(runs){ return new Paragraph({ spacing:{after:120,line:276}, alignment:AlignmentType.JUSTIFIED, children:runs }); }
const body=(t)=>new Paragraph({children:[new TextRun({text:t,font:FONT,size:22})],spacing:{after:120,line:276},alignment:AlignmentType.JUSTIFIED});
const H1=(t)=>new Paragraph({heading:HeadingLevel.HEADING_1,spacing:{before:220,after:120},children:[new TextRun({text:t,bold:true,size:28,font:FONT})]});
const H2=(t)=>new Paragraph({heading:HeadingLevel.HEADING_2,spacing:{before:160,after:100},children:[new TextRun({text:t,bold:true,size:24,font:FONT})]});
const img=(p,w,h,ty)=>new Paragraph({alignment:AlignmentType.CENTER,spacing:{before:120,after:60},children:[new ImageRun({type:ty,data:fs.readFileSync(p),transformation:{width:w,height:h}})]});
const legend=(tag,rest)=>new Paragraph({spacing:{after:170,line:276},children:[new TextRun({text:tag+" ",bold:true,font:FONT,size:20}),new TextRun({text:rest,font:FONT,size:20})]});
const tblCap=(t)=>new Paragraph({spacing:{before:60,after:180,line:252},children:[new TextRun({text:t.split("|")[0],bold:true,font:FONT,size:18}),new TextRun({text:t.split("|").slice(1).join("|"),font:FONT,size:18})]});
const pb=()=>new Paragraph({children:[new PageBreak()]});
function table(rows,widths,fs2){fs2=fs2||16;const total=widths.reduce((a,b)=>a+b,0);
  return new Table({columnWidths:widths,width:{size:total,type:WidthType.DXA},rows:rows.map(function(r,ri){
    return new TableRow({tableHeader:ri===0,children:r.map(function(c,ci){return new TableCell({width:{size:widths[ci],type:WidthType.DXA},
      shading:ri===0?{type:ShadingType.CLEAR,fill:"D9E2F3"}:undefined,
      children:[new Paragraph({spacing:{after:10,line:230},children:[new TextRun({text:String(c),bold:ri===0,font:FONT,size:fs2})]})]});})});}) });}
const J=(f)=>JSON.parse(fs.readFileSync('/home/claude/'+f,'utf8'));
const rez=J('s_rez.json'), red=J('s_red.json'), two=J('s_two.json'), hits=J('s_hits.json'), hits2=J('s_hits2.json'), mets=J('s_metrics.json');
const M='/mnt/user-data/outputs/', MED='/home/claude/ms_unpack/word/media/';
const ch=[]; const push=function(){for(var i=0;i<arguments.length;i++)ch.push(arguments[i]);};

push(new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:80},children:[new TextRun({text:"Supplementary Information (Additional file)",bold:true,size:30,font:FONT})]}));
push(new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:160},children:[new TextRun({text:"Descriptor Selection for Small-Sample QSAR: Selection Bias, Nested Validation and the n:p Ratio Problem , Werle, Rondevaldova, Werle",italics:true,size:20,font:FONT})]}));

push(H1("Contents and data files"));
push(body("This file contains supplementary figures S1 to S6 and supplementary tables S1 to S12. Full descriptor values, correlation matrices and the complete per-strain VIP and effect-size tables are provided as machine-readable spreadsheets (Additional files) rather than reproduced in full here."));
push(body("Additional file 1 - Route B descriptor data (20260706_31MT_LitDeskriptoren_Analyse.xlsx): descriptor definitions and SMARTS, descriptor values for all 31 compounds, Pearson and Spearman matrices, the correlation-based reduction from 25 to 13, and the per-strain PLS-DA and VIP results."));
push(body("Additional file 2 - Route A descriptor matrix (20260706_31MT_5Bac_AllDes_cleaned_CS1-4.xlsx): the cleaned pool of 2,208 descriptors for all 31 compounds, with one sheet per strain documenting the effect-size reduction."));
push(body("Additional file 3 - VIP and effect-size tables (20260710_VIP_EffectSize_Tabellen.xlsx): per-strain VIP and effect-size values, separated by descriptor source, that is the ChemDes Route A1 top-15 set and the literature-derived set."));
push(body("Additional file 4 - Selection-criterion comparison (20260711_Selektion_EffSize_vs_VIP.xlsx): effect-size versus VIP selection at equal model size, containing the underlying values of main-text Table 3 and Table 4."));
push(body("Additional file 5 - Classification metrics (20260720_Klassifikationsmetriken.xlsx): confusion-matrix counts, sensitivity, specificity, balanced accuracy and MCC for every model reported."));
push(body("Additional file 6 - Analysis code and permutation outputs (archive): the Python and RDKit scripts together with the raw permutation results, reproducing every number in this work."));
push(body("Additional file 7 - Source data for the supplementary figures (20260802_FigS5_S6_Rohdaten.xlsx): the LV1 scores underlying Fig. S4 and the two-component scores, loadings and autoscaled descriptor matrices underlying Fig. S5, so that both figures can be redrawn in other software such as MetaboAnalyst."));
push(body("Additional file 8 - Gaucher re-analysis (20260802_Gaucher_Reanalyse.xlsx): full results of applying Routes A1, A1u and A2 to the Smit 2007 data set, together with the 15 highest-ranked m/z values."));

push(pb(),H1("Supplementary figures"));
push(H2("Fig. S1: Single-strain nested permutation example"));
push(img(M+'Fig_nested_permutation_H.I..png',560,250,'png'));
push(legend("Fig. S1.","Detailed fully nested permutation test for Haemophilus influenzae (single-strain example of main-text Fig. 2b). Null distributions of Q2 (left) and AUC (right) from 2,000 permutations of the entire selection-plus-modelling pipeline; the observed value (solid line) lies within the null (dashed: 95th percentile), i.e. not significant."));

push(pb(),H2("Fig. S2: VIP per descriptor"));
push(body("Variable importance in projection (VIP) quantifies how much each descriptor contributes to the latent variable of a fitted PLS-DA model. Apparent VIP is taken from the model fitted on all 31 compounds; nested VIP is the mean across the 31 leave-one-out folds, in which the entire descriptor selection is repeated. The dashed line marks VIP = 1. Because the mean squared VIP equals 1 by construction, descriptors above this line occur in every model, including one fitted to noise, so the line is a ranking aid and not a criterion of validity."));
push(img(M+'Fig_VIP_RouteA_top15.png',430,621,'png'));
push(legend("Fig. S2a.","Route A1, top-15 ChemDes descriptors, all five strains. For each strain the 15 descriptors with the largest effect size after unsupervised correlation grouping of the 2,208 cleaned descriptors are shown, ranked by apparent VIP. The descriptor axis differs between strains because the selection is strain-specific. Apparent and mean nested VIP agree closely, which shows that the VIP ranking is stable across folds; this must not be confused with stability of the selection itself, since which descriptors enter the model varies considerably more than their order once they are in it."));
push(pb(),img(M+'Fig_VIP_RouteB.png',430,621,'png'));
push(legend("Fig. S2b.","Route B, the 13 literature-derived descriptors retained after unsupervised correlation reduction, all five strains. Here the descriptor set is fixed a priori and identical for every strain, so the panels are directly comparable and differences between strains reflect the endpoints rather than the selection. The apparent and nested VIP values coincide because no supervised selection precedes the model."));
push(pb(),img(M+'Fig_VIP_TwoStage.png',430,621,'png'));
push(legend("Fig. S2c.","VIP values for the complete a priori set of 25 literature-derived descriptors, that is Route B before the correlation reduction, all five strains. This is the descriptor set whose performance is compared with the reduced 13-descriptor set in main-text Table 2. Reading this panel against Fig. S2b shows which descriptors the reduction removed and how the importance of the retained ones shifts once their correlated partners are gone. For S. pyogenes, the one endpoint that the reduction demonstrably harms, several descriptors carrying intermediate VIP here are absent from the reduced set."));

push(pb(),H2("Fig. S3: VIP curves"));
push(body("In all three panels the descriptors are ordered along the x-axis by their mean VIP across the five strains, from highest to lowest. This ordering is deliberately strain-independent, so that the curve of an individual strain can be read as a deviation from the common profile rather than from its own ranking; a strain whose curve falls steeply from left to right shares the common ordering, whereas a strain whose curve is flat or irregular weights the descriptors differently from the others."));
push(img(M+'Fig_VIPcurve_RouteB.png',560,255,'png'));
push(legend("Fig. S3a.","VIP curve for the 13-descriptor Route B set, all five strains overlaid on a common descriptor axis sorted by mean VIP. Because the set is fixed a priori and identical for every strain, the curves are directly comparable. Two patterns stand out. S. aureus and S. pneumoniae run almost in parallel, so the same descriptors drive both endpoints, which is consistent with these being the two endpoints that remain significant under every route. S. pyogenes and H. influenzae depart clearly from that pattern and from each other, S. pyogenes weighting the polarity and ring terms most heavily. The descriptors that matter therefore differ between bacteria, which is a biological statement and not an artefact of selection, since no selection took place."));
push(img(M+'Fig_VIPcurve_TwoStage.png',560,255,'png'));
push(legend("Fig. S3b.","VIP curve for the complete a priori set of 25 literature-derived descriptors, all strains overlaid and sorted by mean VIP. The grouping seen in Fig. S3a is already present before the correlation reduction, which shows that it is a property of the endpoints and not a consequence of which descriptors happened to survive that reduction."));
push(pb(),img(M+'Fig_VIPcurve_RouteA1_top15.png',600,273,'png'));
push(legend("Fig. S3c.","VIP curve for the Route A1 top-15 descriptors, drawn in the same way as Fig. S3a and Fig. S3b. The descriptor axis is the union of the five strain-specific top-15 sets, so a gap means that the descriptor was not selected for that strain. The many gaps are themselves the result: unlike the fixed literature set, the data-driven selection retains largely different descriptors for different endpoints, which is the instability discussed in the main text."));

push(pb(),H2("Fig. S4: One-component (LV1) apparent class separation"));
push(img(M+'Fig_LV1_boxjitter_all.png',620,177,'png'));
push(legend("Fig. S4.","One-component (LV1) apparent class separation for all five strains, Route A1 top-15, shown as a box plot with the individual molecules overlaid. The box gives the median and interquartile range and every compound is plotted, which matters at these very small class sizes. The LV1 axis is oriented so that active compounds score higher. The separation shown here is apparent, that is obtained from a model fitted to all compounds, and should be read together with the nested validation in Fig. 3b."));

push(pb(),H2("Fig. S5: Two-dimensional PLS-DA score plots"));
push(img(M+'Fig_PLSDA_scoreplots_all.png',600,380,'png'));
push(legend("Fig. S5.","MetaboAnalyst-style two-component PLS-DA score plots (samples) for all five strains, Route A1 top-15: active (red) and non-active (green) with 95% confidence ellipses. The apparent separation is visually convincing for every strain, including those that fail nested validation, and is therefore not evidence of predictive ability."));

push(pb(),H2("Fig. S6: Full null distributions of the univariate hit count, all three reduction schemes"));
push(img(M+'Fig_correlation_inflation.png',575,380,'png'));
push(legend("Fig. S6.","How the dependence between the univariate tests changes with the reduction scheme, shown as the full null distribution rather than the summary inflation factor of main-text Fig. 5a. For each strain the grey histogram is the null distribution of the number of descriptors passing the effect-size filter over 1,000 label permutations; the solid red line is the observed count with the true labels, the dotted line the naive expectation of 5% of the tests performed, and the dashed line the 95th percentile of the null. (a) Route A2, in which the filter is applied to all 2,208 descriptors: the null distribution is extremely broad and the observed counts fall inside it. (b) Route A1, in which the 2,208 descriptors are first grouped to 182 representatives chosen by effect size. (c) Route A1u, in which the representative is chosen by centrality so that the grouping is label-free: the null distribution is now narrow and close to the naive expectation, and the excess observed for S. aureus and S. pneumoniae becomes interpretable (compare main-text Fig. 5b, which shows the Route A1u panel as a boxplot)."));

push(pb(),H1("Supplementary tables"));
push(H2("Table S1: Route B descriptor definitions and conventions"));
var s1=[["Descriptor","Class","Method / RDKit","SMARTS","Convention / note"]].concat(rez.map(function(r){return [r["Deskriptor"],r["Gruppe"],r["Methode/RDKit"],r["SMARTS"]||"not applicable",r["Konvention/Notiz"]];}));
push(table(s1,[1500,1150,1650,2050,2450],13));
push(tblCap("Table S1|  The 25 a priori literature-derived descriptors, computed with RDKit from neutral canonical SMILES (no 3D optimisation, no pH correction)."));

push(H2("Table S2: Unsupervised correlation-based reduction (25 to 13)"));
var s2=[["Representative (kept)","Type","Redundant partners (|r| > 0.7)","Cluster size"]].concat(red.map(function(r){return [r["Repräsentant (behalten)"],r["Typ"],r["redundante Partner (|r|>0.7)"]||"none (singleton)",r["Clustergröße"]];}));
push(table(s2,[1900,1700,3600,900],15));
push(tblCap("Table S2|  Correlation-based feature grouping of the 25 Route-B descriptors (Pearson or Spearman |r| > 0.7). Rare binary indicators (< 4 positives) kept as singletons. The 13 representatives form the modelling set."));

push(pb(),H2("Table S3: Complete results across the three analyses"));
push(H2("Table S3a: Route A1 (grouping first), full precision"));
var s3a=[
 ["Strain","app Q2","app AUC","nested Q2","nested AUC","p(Q2)","p(AUC)"],
 ["H. influenzae","0.40","0.90","+0.02","0.70","0.080","0.077"],
 ["S. aureus","0.51","0.97","+0.10","0.81","0.032","0.032"],
 ["S. pneumoniae","0.56","0.97","+0.30","0.89","0.003","0.002"],
 ["S. pyogenes","0.49","0.94","-0.13","0.65","0.205","0.135"],
 ["P. aeruginosa †","0.66","1.00","+0.41","0.99","artefact","artefact"]];
push(table(s3a,[2200,1050,1050,1200,1200,1000,1000],12));
push(tblCap("Table S3a|  Route A1 (unsupervised correlation grouping of all 2,208 descriptors, then the effect-size filter): apparent Q2/AUC from selection performed once on all compounds; nested Q2/AUC from selection repeated inside every leave-one-out fold with whole-pipeline permutation (2,000 permutations); permutation p-values for both metrics. One latent variable throughout. † P. aeruginosa (3 actives) is an artefact."));

push(H2("Table S3b: Route A2 (effect-size filter first), full precision"));
var s3b=[
 ["Strain","app Q2","app AUC","nested Q2","nested AUC","p(Q2)","p(AUC)"],
 ["H. influenzae","0.46","0.94","-0.14","0.70","0.154","0.080"],
 ["S. aureus","0.52","0.97","+0.15","0.82","0.013","0.023"],
 ["S. pneumoniae","0.55","0.97","+0.19","0.87","0.009","0.006"],
 ["S. pyogenes","0.53","0.95","-0.11","0.71","0.135","0.056"],
 ["P. aeruginosa †","0.66","1.00","+0.40","0.99","artefact","artefact"]];
push(table(s3b,[2200,1050,1050,1200,1200,1000,1000],12));
push(tblCap("Table S3b|  Route A2 (effect-size filter applied to all 2,208 descriptors first, then correlation grouping of the survivors), same layout and same 2,000 permutations as Table S3a. † P. aeruginosa (3 actives) is an artefact."));

push(H2("Table S3c: Route B, full precision"));
var s3c=[
 ["Strain","Q2 (nLV=1)","AUC (nLV=1)","p(Q2)","p(AUC)"],
 ["H. influenzae","0.29","0.85","0.002","0.001"],
 ["S. aureus","0.39","0.82","<0.001","0.016"],
 ["S. pneumoniae","0.35","0.85","0.001","0.008"],
 ["S. pyogenes","0.06","0.62","0.034","0.127"],
 ["P. aeruginosa †","0.82","1.00","artefact","artefact"]];
push(table(s3c,[2200,1450,1450,1050,1050],12));
push(tblCap("Table S3c|  Route B (13 literature-derived descriptors, unsupervised reduction, one latent variable). Because no supervised selection precedes the model, apparent and nested values coincide and a simple (non-nested) permutation test is already honest. † P. aeruginosa (3 actives) is an artefact."));

push(H2("Table S4: An additional VIP variant applied to the 25 literature descriptors"));
var s4=[["Strain","nLV","Q2 naive","AUC naive","Q2 nested","AUC nested","# VIP > 1"]].concat(
  two.map(function(r){return [r["Erreger"],r["nLV"],r["Q2_naiv"],r["AUC_naiv"],r["Q2_genestet"],r["AUC_genestet"],r["nVIP_gt1_gesamt"]];}));
push(table(s4,[1500,700,1250,1200,1250,1200,1000],15));
push(tblCap("Table S4|  A further two-stage VIP variant, retained here for completeness. It differs from the procedure reported in the main text (Fig. 1c, Table 4) in two respects: it starts from the 25 literature-derived descriptors rather than the Route A1 set, and it selects by the conventional threshold of VIP above 1 rather than by a fixed model size, with the number of latent variables in the first model chosen to maximise Q2. Both of these choices are argued against in the main text, the first because it mixes the two routes and the second because the VIP above 1 threshold depends on how many descriptors entered the model (Xu and colleagues, reference 4 of the main text); the variant is therefore not part of the main analysis. It is shown because it reproduces the same qualitative result independently: the naive Q2, obtained with the VIP selection performed once on all compounds, is consistently higher than the nested Q2, and that difference is the optimistic bias introduced by label-dependent selection outside cross-validation."));

push(H2("Table S5: Chance descriptor counts under the effect-size filter for the three reduction schemes"));
var s5=[["Strain","A2 obs","A2 null mean","A2 null SD","A2 p","A1 obs","A1 null mean","A1 null SD","A1 p","A1u obs","A1u null mean","A1u null SD","A1u p"]].concat(
  hits2.map(function(r){return [r.strain,r.a2o,r.a2m,r.a2s,r.a2p,r.a1o,r.a1m,r.a1s,r.a1p,r.auo,r.aum,r.aus,r.aup];}));
push(table(s5,[1150,600,750,650,550,600,750,650,550,650,780,680,580],11));
push(tblCap("Table S5|  Number of descriptors passing the effect-size filter for the true labels (obs) and under 1,000 label permutations (null mean, null standard deviation, and the permutation p-value for the observed count). A2 applies the filter to all 2,208 descriptors; A1 and A1u apply it to the 182 representatives obtained by unsupervised correlation grouping, with the representative chosen by effect size (A1) or by centrality (A1u). Under independence the standard deviation would be about 10.2 for 2,208 tests and about 2.9 for 182 tests, so the observed dependence inflates it about ten-fold under A2, about 3.3-fold under A1 and about 1.7-fold under A1u. Only under A1u is the observed count interpretable, and it is then significantly above chance for S. aureus and S. pneumoniae (main-text Fig. 5)."));

push(pb(),H2("Table S6: Classification metrics for all reported models"));
var s6=[["Strain","Model","Desc.","TP","FN","FP","TN","Sens.","Spec.","Bal. acc.","MCC"]].concat(
  mets.map(function(r){return [r.Strain,r.Model,r.Descriptors,r.TP,r.FN,r.FP,r.TN,r.Sensitivity,r.Specificity,r.BalancedAccuracy,r.MCC];}));
push(table(s6,[1250,1500,600,450,450,450,450,650,650,750,650],11));
push(tblCap("Table S6|  Confusion-matrix counts and derived classification metrics for every model reported in this study, obtained from the leave-one-out predictions at a decision threshold of 0.5. Apparent denotes selection performed once on all compounds; nested denotes selection repeated inside every fold. A1 groups the descriptors before filtering and chooses the group representative by effect size, A1u is the same with the representative chosen by centrality, and A2 applies the effect-size filter first. B is the literature-derived set. Sensitivity is the fraction of active compounds recovered, specificity the fraction of inactive compounds correctly classified, balanced accuracy their mean, and MCC the Matthews correlation coefficient, which ranges from -1 to +1 with 0 corresponding to chance. The drop from the apparent to the nested variant is visible in the Matthews coefficient just as it is in Q2."));

push(pb(),H2("Table S7: The same three routes applied to an independent published data set, with the original results for comparison"));
var s8=[
 ["", "A1, 1 LV", "A1u, 1 LV", "A2, 1 LV", "A1, 2 LV", "A1u, 2 LV", "A2, 2 LV", "Smit 2007"],
 ["Classifier", "PLS-DA", "PLS-DA", "PLS-DA", "PLS-DA", "PLS-DA", "PLS-DA", "PCDA"],
 ["Variables entering the model", "68", "60", "69", "68", "60", "69", "590"],
 ["Latent variables / components", "1", "1", "1", "2", "2", "2", "15 (7-20)"],
 ["Apparent Q2", "0.71", "0.76", "0.74", "0.80", "0.81", "0.80", "n.r."],
 ["Apparent AUC", "0.99", "0.99", "0.99", "1.00", "0.99", "0.99", "n.r."],
 ["Nested Q2 (p)", "0.39 (<0.001)", "0.45 (<0.001)", "0.38 (<0.001)", "0.50 (0.001)", "0.58 (0.001)", "0.45 (0.001)", "n.r."],
 ["Nested AUC (p)", "0.87 (<0.001)", "0.89 (<0.001)", "0.88 (<0.001)", "0.93 (0.001)", "0.96 (0.001)", "0.92 (0.001)", "n.r."],
 ["Nested sensitivity", "0.79", "0.74", "0.79", "0.84", "0.89", "0.84", "0.89"],
 ["Nested specificity", "0.85", "0.85", "0.85", "0.85", "0.90", "0.85", "0.90"],
 ["Nested balanced accuracy", "0.82", "0.79", "0.82", "0.85", "0.90", "0.85", "0.90"],
 ["Nested MCC", "0.64", "0.59", "0.64", "0.69", "0.79", "0.69", "0.79"],
 ["Misclassified (of 39)", "7", "8", "7", "6", "4", "6", "4"],
];
push(table(s8,[1900,930,930,930,930,930,930,980],9));
push(tblCap("Table S7|  Routes A1, A1u and A2 applied unchanged to the SELDI-TOF serum proteomics data of Smit and colleagues (2007): 20 healthy controls, 19 patients with Gaucher disease and 590 m/z variables per sample. The pipeline, the leave-one-out cross-validation and the 2,000 whole-pipeline permutations are exactly as for the monoterpenoid data; only the input matrix was exchanged, and autoscaling was performed inside every fold on the training samples only. The first three columns use a single latent variable, the conservative setting used throughout this work. The next three use two latent variables, which is the number that maximises the cross-validated Q2 for every route; the component number is fixed a priori in these columns, so the permutation test remains valid and every value shown is significant at p < 0.002. Choosing the component number on the data costs remarkably little here: repeating the analysis with the number selected in an inner cross-validation loop, that is a doubly nested design of the kind used in the original publication, gives Q2 = 0.48 for A1, 0.57 for A1u and 0.44 for A2, within 0.02 of the fixed two-component values, and the inner loop selects two components in the majority of folds. The last column reproduces the figures given by Smit and colleagues for their double cross-validated principal component discriminant analysis, whose component number was tuned in an inner loop and ranged from 7 to 20 with 15 selected on the full data. That work reported neither Q2 nor AUC, which is marked n.r.; the balanced accuracy and the Matthews coefficient shown for it are derived here from the confusion matrix implied by their reported 4 misclassifications, namely 2 controls and 2 patients. Route A1u with two components reproduces the published sensitivity, specificity and misclassification count exactly, using a single tuned parameter instead of up to twenty components."));

push(pb(),H2("Additional analysis: does a second latent variable help?"));
push(rich([T("The single latent variable used throughout was chosen a priori as the conservative option, and the Gaucher comparison raises the fair question of whether it is too conservative here as well. We therefore repeated the analysis with two and three components and, to separate a genuine gain from tuning optimism, also with the component number selected by an inner cross-validation on the training compounds only (Table S8). The two routes answer the question differently, and the difference is instructive."),
 T(" Under Route A1 a second component is harmful. The nested Q\u00b2 falls for every evaluable endpoint except "), I("S. pyogenes"), T(", and for "), I("H. influenzae"),
 T(" it collapses from 0.02 to -0.43. The inner loop, given a free choice between one and four components, selects a single component in the majority of folds for all four evaluable endpoints, and the doubly nested estimate never exceeds the fixed one-component value. The a priori choice is therefore confirmed by an honest tuner rather than merely assumed.")]));
push(rich([T("Under Route B the opposite holds. With the fixed 13-descriptor set, three components raise the honest Q\u00b2 for "),
 I("H. influenzae"), T(" from 0.29 to 0.46 (p = 0.002) and for "), I("S. pneumoniae"), T(" from 0.35 to 0.44 (p = 0.002), while "),
 I("S. aureus"), T(" gains slightly and "), I("S. pyogenes"),
 T(" loses. The doubly nested estimates confirm that these gains survive an honest choice of the component number. A plausible reading is that the supervised selection of Route A1 has already condensed what signal there is into the first component, so that further components can only fit noise, whereas Route B, which spends no degrees of freedom on selection, still has structure left for a second and third component to extract. Both settings are therefore reported for Route B in Table 1a and Table 1b: the one-component model for direct comparability with Route A, and the three-component model because it is the better description of what the interpretable descriptor set can achieve. For Route A only the one-component model is reported, since that is the setting an honest tuner selects there.")]));
push(rich([T("One further observation belongs here. The nested Q\u00b2 of "), I("P. aeruginosa"),
 T(" rises monotonically with model capacity, from 0.41 with one component through 0.59 with two to 0.73 when the component number is chosen freely. A genuine signal does not usually improve without limit as capacity is added; an artefact produced by three active compounds does exactly that, because each additional component can be spent on fitting them. This behaviour is an independent indication that the endpoint is not evaluable, alongside the argument from its class size.")]));
push(H2("Table S8: Effect of the number of latent variables in the monoterpenoid data"));
var s9=[
 ["", "H. influenzae", "S. aureus", "S. pneumoniae", "S. pyogenes", "P. aeruginosa +"],
 ["Route A1, nested Q2, 1 LV", "0.02", "0.10", "0.30", "-0.13", "0.41"],
 ["Route A1, nested Q2, 2 LV", "-0.43", "-0.01", "0.17", "-0.04", "0.59"],
 ["Route A1, nested Q2, 3 LV", "-0.33", "-0.20", "0.22", "-0.17", "0.72"],
 ["Route A1, doubly nested Q2", "-0.07", "0.06", "0.30", "-0.13", "0.73"],
 ["Route A1, components chosen inside", "1 (1-4)", "1 (1-3)", "1 (1-1)", "1 (1-1)", "4 (2-4)"],
 ["Route B, Q2, 1 LV", "0.29", "0.39", "0.35", "0.06", "0.82"],
 ["Route B, Q2, 2 LV", "0.35", "0.37", "0.36", "-0.01", "0.86"],
 ["Route B, Q2, 3 LV (p)", "0.46 (0.002)", "0.42 (<0.001)", "0.44 (0.002)", "-0.01 (0.048)", "0.96 (0.002)"],
 ["Route B, doubly nested Q2", "0.44", "0.31", "0.37", "-0.00", "0.89"],
 ["Route B, components chosen inside", "4 (3-4)", "4 (1-4)", "3 (1-4)", "1 (1-4)", "4 (2-4)"],
];
push(table(s9,[2500,1400,1250,1400,1300,1400],10));
push(tblCap("Table S8|  Nested Q2 as a function of the number of latent variables, for the data-driven Route A1 and the a priori Route B. For Route A1 the descriptor selection is repeated inside every fold as everywhere in this work; for Route B the descriptor set is fixed, so no nesting of a selection step is required. The rows labelled doubly nested give the estimate obtained when the number of components is chosen by an inner cross-validation on the training compounds only, and the following row reports the median and range of the numbers that inner loop selected across the 31 outer folds. Permutation p-values from 1,000 permutations are given for the Route B three-component models. The two routes behave in opposite ways. Under Route A1 a second component is harmful for every evaluable endpoint except S. pyogenes, dramatically so for H. influenzae, and the inner loop selects a single component in the majority of folds, so the a priori choice of one component made throughout this work is confirmed by an honest tuner rather than merely assumed. Under Route B, where no supervised selection precedes the model, additional components help three of the four evaluable endpoints, raising Q2 for H. influenzae from 0.29 to 0.46 and for S. pneumoniae from 0.35 to 0.44, and the doubly nested values confirm that this gain is not tuning optimism. The dagger marks P. aeruginosa, whose apparent advantage grows monotonically with model capacity, from 0.41 through 0.59 to 0.73, which is the behaviour expected of an artefact rather than of a signal."));

push(pb(),H2("Additional analysis: the cross-validated c-statistic is itself conservative"));
push(body("One property of our reporting deserves explicit mention. The area under the ROC curve is identical to the concordance statistic, or c-statistic, of the prognostic modelling literature, and the leave-one-out cross-validated c-statistic is known to be biased towards zero, most strongly for estimators that shrink predictions towards the observed event fraction and for small samples with rare events [39]. Our setting has all three features. Recomputing the c-statistic by leave-pair-out, in which one active and one inactive compound are removed together and the entire pipeline is refitted, shows that the bias is present and uniformly conservative: leave-one-out understates it by 0.038 on average for Route A1 and by up to 0.090 for individual endpoints (Table S9). Every AUC reported here is therefore, if anything, too pessimistic, which strengthens rather than weakens the significant results and softens the negative ones. Three considerations nevertheless argue against changing the reported scheme. The permutation test compares the observed statistic with a null distribution generated by identical machinery, so both are affected equally and the test remains valid. Q2, our primary metric, is algebraically one minus the ratio of the model Brier score to that of the null model, and the leave-one-out Brier score is nearly unbiased, so Q2 is largely unaffected. And leave-one-out is the only scheme for which the whole-pipeline permutation test is computationally feasible at all, as set out in the legend of Table S9. The one caveat we would flag is that a bias of this kind is not identical across routes, so comparisons of AUC between routes carry a small distortion; in our data it runs against Route B, whose advantage is therefore understated."));
push(H2("Table S9: Resampling scheme and the bias of the cross-validated c-statistic"));
var s10=[
 ["","H. influenzae","S. aureus","S. pneumoniae","S. pyogenes","P. aeruginosa +"],
 ["Route A1, LOO AUC", "0.70", "0.81", "0.89", "0.65", "0.99"],
 ["Route A1, leave-pair-out c-statistic", "0.75", "0.85", "0.93", "0.71", "0.99"],
 ["Route A1, difference", "+0.05", "+0.04", "+0.04", "+0.06", "+0.00"],
 ["Route A1u, LOO AUC", "0.55", "0.78", "0.81", "0.60", "0.87"],
 ["Route A1u, leave-pair-out c-statistic", "0.64", "0.83", "0.87", "0.66", "0.94"],
 ["Route A1u, difference", "+0.09", "+0.05", "+0.05", "+0.06", "+0.07"],
 ["Route A2, LOO AUC", "0.70", "0.82", "0.87", "0.71", "0.99"],
 ["Route A2, leave-pair-out c-statistic", "0.74", "0.88", "0.91", "0.78", "0.99"],
 ["Route A2, difference", "+0.04", "+0.05", "+0.04", "+0.07", "+0.00"],
 ["Route B, LOO AUC", "0.85", "0.82", "0.85", "0.62", "1.00"],
 ["Route B, leave-pair-out c-statistic", "0.88", "0.90", "0.91", "0.71", "1.00"],
 ["Route B, difference", "+0.02", "+0.08", "+0.06", "+0.09", "+0.00"],
 ["Route A1, leave-2-out Q2", "0.04", "0.11", "0.29", "-0.04", "0.38"],
 ["Route A1, leave-3-out Q2", "0.04", "0.10", "0.29", "-0.01", "0.36"],
 ["Route A1, leave-one-out Q2", "0.02", "0.10", "0.30", "-0.13", "0.41"],
];
push(table(s10,[2750,1250,1150,1300,1200,1300],10));
push(tblCap("Table S9|  Effect of the resampling scheme on the reported performance. The upper part compares the leave-one-out AUC with the leave-pair-out c-statistic for all four routes. Leave-pair-out removes one active and one inactive compound together, refits the entire pipeline on the remaining 29 compounds and scores the pair as concordant if the active receives the higher prediction; it was proposed because the leave-one-out c-statistic is biased towards zero, most strongly for shrinkage estimators and small, imbalanced samples [39]. The bias is present in our data and is uniformly conservative: leave-one-out understates the c-statistic by 0.038 on average for Route A1, 0.041 for Route A2, 0.051 for Route B and 0.065 for Route A1u, with a maximum of 0.090. Every reported AUC is therefore, if anything, too pessimistic. The lower part gives the nested Q2 of Route A1 under leave-one-out, leave-2-out and leave-3-out. The differences are at most 0.03 for four of the five endpoints, so the conclusions do not depend on the resampling scheme; the exception is S. pyogenes, which improves from -0.13 to -0.01 without approaching significance. Note that P. aeruginosa moves in the opposite direction, its apparent advantage falling as more compounds are held out, and that 1.9 % of its leave-3-out folds cannot be modelled at all because too few actives remain. Leave-one-out is nevertheless retained throughout this work for a practical reason: only for it is the whole-pipeline permutation test computationally feasible. Repeating the 2,000 permutations for five strains would require about 0.3 million model fits under leave-one-out, about 4.7 million under leave-2-out and about 45 million under leave-3-out, the last corresponding to roughly 470 hours at the speed measured here. The values in the lower part of this table are therefore point estimates without permutation p-values. The dagger marks P. aeruginosa as an artefact."));

push(pb(),H2("Table S10: How much signal survives at our sample sizes"));
var s11=[
 ["", "19 cases", "11 cases", "10 cases", "6 cases", "5 cases", "4 cases"],
 ["Matched monoterpenoid endpoint", "none (full data)", "S. pyogenes (11 / 31)", "H. influenzae (10 / 31)", "S. pneumoniae (6 / 31)", "S. aureus (5 / 31), count", "S. aureus (5 / 31), fraction"],
 ["Cases drawn (of 19)", "19", "11", "10", "6", "5", "4"],
 ["Controls retained (of 20)", "20", "20", "20", "20", "20", "20"],
 ["Subsample size n", "39", "31", "30", "26", "25", "24"],
 ["Active fraction, subsample", "48.7 %", "35.5 %", "33.3 %", "23.1 %", "20.0 %", "16.7 %"],
 ["Active fraction, our endpoint", "-", "35.5 %", "32.3 %", "19.4 %", "16.1 %", "16.1 %"],
 ["Random draws", "1", "150", "150", "150", "150", "150"],
 ["Nested Q2, median", "0.39", "0.37", "0.34", "0.21", "0.15", "0.07"],
 ["Nested Q2, 5th to 95th percentile", "single value", "0.16 to 0.50", "0.09 to 0.51", "-0.14 to 0.44", "-0.30 to 0.40", "-0.42 to 0.38"],
 ["Nested AUC, median", "0.87", "0.86", "0.85", "0.81", "0.79", "0.76"],
 ["Draws with Q2 above 0", "100 %", "100 %", "98 %", "81 %", "75 %", "65 %"],
 ["Draws with Q2 above 0.30", "100 %", "66 %", "59 %", "27 %", "17 %", "12 %"],
];
push(table(s11,[2450,1150,1000,1000,1000,1050,1050],9));
push(tblCap("Table S10|  Detection limit of the analysis pipeline, established on data whose signal is independently confirmed. Procedure: from the 39 Gaucher samples of Smit and colleagues, k of the 19 patients were drawn at random without replacement while all 20 controls were retained unchanged, giving a subsample of n = k + 20; Route A1 with one latent variable was then rerun on that subsample exactly as in the main analysis, with the descriptor selection repeated inside every leave-one-out fold. This was repeated for 150 independent draws per column, and the table reports the median, the 5th to 95th percentile and the proportion of draws exceeding a given value. The column with 19 cases is the full data set and therefore a single deterministic value. The threshold 0.30 is approximately the best nested Q2 obtained anywhere in the monoterpenoid data. Because only 20 controls exist, a subsample cannot reproduce both the number of actives and the sample size of our endpoints at once. Drawing controls with replacement would remove that constraint but must not be done here: with 26 controls drawn from 20, about 11 are duplicates on average, and whenever leave-one-out holds out a duplicated control its identical twin remains in the training set, so the model has already seen the held-out sample. We verified the size of this effect rather than assuming it: bootstrapping the controls to n = 31 raises the median nested Q2 by 0.10 for the S. aureus configuration and by 0.13 for S. pneumoniae, and lifts the proportion of draws with a positive Q2 from 75 % to 90 %. That is leakage of exactly the kind this study is about, acting in the direction that would flatter our argument, and the subsamples were therefore drawn without replacement throughout. The columns with 11, 10, 6 and 5 cases match the number of actives, and the last column instead matches the active fraction of S. aureus, which the count-matched column overstates (20.0 % against 16.1 %). Where the two differ, the fraction-matched column is the more demanding and the more realistic comparison: its median Q2 of 0.07 is less than half the count-matched value. The result cuts both ways and is reported here for that reason. With ten or eleven cases a strong signal remains visible in essentially every draw, so the negative findings for H. influenzae and S. pyogenes cannot be attributed to a lack of power. With four to six cases the same signal yields a negative Q2 in a fifth to a third of the draws and exceeds 0.30 in at most a quarter, so the positive findings for S. aureus and S. pneumoniae rest on estimates of very large sampling variability. Two further differences from our own data should be kept in mind: the subsamples retain 590 variables against our 2,208, giving a more favourable sample-to-variable ratio, and the underlying effect is a diagnosed disease state, almost certainly stronger than any structure-activity relationship in a congeneric series. These figures are therefore an optimistic bound on what is detectable at these class sizes."));

push(pb(),H2("Table S11: Why the gap between R2 and Q2 does not detect selection leakage"));
var s12=[
 ["", "H. influenzae", "S. aureus", "S. pneumoniae", "S. pyogenes", "P. aeruginosa +"],
 ["Route A1, R2 (fit on all compounds)", "0.53", "0.63", "0.65", "0.60", "0.76"],
 ["Route A1, Q2 apparent", "0.40", "0.51", "0.56", "0.49", "0.66"],
 ["Route A1, R2 minus Q2 apparent", "0.13", "0.12", "0.09", "0.11", "0.10"],
 ["Route A1, Q2 nested", "0.02", "0.10", "0.30", "-0.13", "0.41"],
 ["Route A1, Q2 apparent minus Q2 nested", "0.38", "0.41", "0.26", "0.62", "0.24"],
 ["Route B, R2 (fit on all compounds)", "0.52", "0.58", "0.55", "0.46", "0.87"],
 ["Route B, Q2", "0.29", "0.39", "0.35", "0.06", "0.82"],
 ["Route B, R2 minus Q2", "0.23", "0.19", "0.20", "0.40", "0.05"],
];
push(table(s12,[2900,1200,1100,1250,1150,1250],10));
push(tblCap("Table S11|  Apparent fit (R2), cross-validated performance (Q2) and the difference between them, for the data-driven Route A1 and the leakage-free Route B, all with one latent variable. R2 is the coefficient of determination of the model fitted to all 31 compounds using the descriptors selected on those same compounds. A large gap between R2 and Q2 is a recognised warning sign of overfitting, and a threshold of about 0.3 is often quoted. By that criterion every Route A1 model here would pass comfortably, with gaps of 0.09 to 0.13, yet the same models lose between 0.26 and 0.62 of their Q2 when the descriptor selection is repeated inside every fold. For S. pyogenes the discrepancy is extreme: a gap of 0.11 accompanies a nested Q2 of -0.13, that is a model that predicts worse than assigning every compound the class mean. The comparison with Route B makes the reason visible. There the descriptors are fixed a priori, no selection precedes the model, and the gaps are larger, 0.19 to 0.40, although these are the honestly predictive models. In other words the criterion is not merely insensitive here but inverted, because R2 and Q2 are both computed on a descriptor set that was already chosen using the same labels; the difference between them reflects only the modelling step, while the leakage resides in the selection step that precedes both. The dagger marks P. aeruginosa as an artefact."));

push(pb(),H2("Table S12: A wrapper method for comparison, and how much more it inflates its own estimate"));
var s13=[
 ["", "H. influenzae", "S. aureus"],
 ["Descriptors offered to the search", "182", "182"],
 ["Descriptors selected (single run on all data)", "26", "24"],
 ["Q2 reported by the search itself", "0.64", "0.72"],
 ["Q2, LOO with that set held fixed", "0.44", "0.70"],
 ["AUC, LOO with that set held fixed", "0.88", "0.97"],
 ["Q2, search repeated inside every fold", "-0.59", "-0.08"],
 ["AUC, search repeated inside every fold", "0.51", "0.68"],
 ["Collapse (self-reported minus honest Q2)", "1.23", "0.80"],
 ["Descriptors per fold, median (range)", "23 (13 to 29)", "23 (16 to 33)"],
 ["For comparison, Route A1 nested Q2", "0.02", "0.10"],
 ["For comparison, Route A1 collapse", "0.42", "0.42"],
];
push(table(s13,[4100,1700,1700],10));
push(tblCap("Table S12|  A genetic algorithm applied to the same data, for comparison with the filter used throughout this work. The implementation follows the design described for the OPERA models: a binary chromosome encoding presence or absence of each descriptor, a population of 30, crossover probability 0.5, mutation probability 0.01, 50 generations, and a fitness function that maximises cross-validated Q2 while penalising model size. The search was offered the 182 label-free correlation-group representatives rather than the full pool, because the latter spans about 10 to the power 665 subsets. The row labelled Q2 reported by the search itself is the fitness of the winning chromosome, that is the number such a procedure would ordinarily publish. The next two rows evaluate that same fixed descriptor set by leave-one-out. The rows below repeat the entire search inside every leave-one-out fold, which is the honest estimate. The collapse is two to three times larger than the 0.42 seen for the effect-size filter, and it is worst for S. aureus, the one endpoint for which our filter finds a reproducible signal. Note also that the search selected a different subset in every fold, of median size 23. This is not a criticism of the cited work, which validated its models on a held-out external test set and is therefore unaffected; it is an illustration that the amount of search freedom, not the choice of algorithm, governs how far a procedure inflates its own estimate. A full replication of that protocol, with 100 runs in each of two rounds, would require about 79 hours nested and roughly 18 years if combined with 2,000 permutations, so the values here come from a single run per fold and are indicative rather than definitive."));

const doc=new Document({creator:"Werle et al.",styles:{default:{document:{run:{font:FONT,size:22}}}},
  sections:[{properties:{page:{margin:{top:1440,bottom:1440,left:1080,right:1080}}},children:ch}]});
Packer.toBuffer(doc).then(function(b){fs.writeFileSync(M+"Werle_JCheminform_SupplementaryInformation_2026-08-10.docx",b);console.log("written");});
