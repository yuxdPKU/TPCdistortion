#ifdef __CLING__
#pragma cling optimize(0)
#endif
void pipi_nopull()
{
//=========Macro generated from canvas: massFitCanvas/massFitCanvas
//=========  (Tue May 27 09:19:27 2025) by ROOT version 6.32.06
   TCanvas *massFitCanvas = new TCanvas("massFitCanvas", "massFitCanvas",0,0,800,800);
   gStyle->SetOptTitle(0);
   massFitCanvas->Range(0,0,1,1);
   massFitCanvas->SetFillColor(0);
   massFitCanvas->SetBorderMode(0);
   massFitCanvas->SetBorderSize(2);
   massFitCanvas->SetTickx(1);
   massFitCanvas->SetTicky(1);
   massFitCanvas->SetLeftMargin(0.16);
   massFitCanvas->SetRightMargin(0.05);
   massFitCanvas->SetTopMargin(0.08);
   massFitCanvas->SetBottomMargin(0.16);
   massFitCanvas->SetFrameBorderMode(0);
   
   TH1D *frame_K_S0_mass_1aa3fce0__1 = new TH1D("frame_K_S0_mass_1aa3fce0__1","A RooPlot of \"mass\"",100,0.4,0.6);
   frame_K_S0_mass_1aa3fce0__1->SetBinContent(1,5482.852);
   frame_K_S0_mass_1aa3fce0__1->SetMaximum(5482.852);
   frame_K_S0_mass_1aa3fce0__1->SetEntries(1);
   frame_K_S0_mass_1aa3fce0__1->SetDirectory(nullptr);
   frame_K_S0_mass_1aa3fce0__1->SetStats(0);
   frame_K_S0_mass_1aa3fce0__1->SetMarkerStyle(4);
   frame_K_S0_mass_1aa3fce0__1->SetMarkerSize(0.02);
   frame_K_S0_mass_1aa3fce0__1->GetXaxis()->SetTitle("m(#pi^{+}#pi^{-}) [GeV]");
   frame_K_S0_mass_1aa3fce0__1->GetXaxis()->SetNdivisions(605);
   frame_K_S0_mass_1aa3fce0__1->GetXaxis()->SetLabelFont(42);
   frame_K_S0_mass_1aa3fce0__1->GetXaxis()->SetLabelSize(0.05);
   frame_K_S0_mass_1aa3fce0__1->GetXaxis()->SetTitleSize(0.05);
   frame_K_S0_mass_1aa3fce0__1->GetXaxis()->SetTitleOffset(1.4);
   frame_K_S0_mass_1aa3fce0__1->GetXaxis()->SetTitleFont(42);
   frame_K_S0_mass_1aa3fce0__1->GetYaxis()->SetTitle("Candidates / (4 MeV)");
   frame_K_S0_mass_1aa3fce0__1->GetYaxis()->SetLabelFont(42);
   frame_K_S0_mass_1aa3fce0__1->GetYaxis()->SetLabelSize(0.05);
   frame_K_S0_mass_1aa3fce0__1->GetYaxis()->SetTitleSize(0.05);
   frame_K_S0_mass_1aa3fce0__1->GetYaxis()->SetTitleOffset(1.4);
   frame_K_S0_mass_1aa3fce0__1->GetYaxis()->SetTitleFont(42);
   frame_K_S0_mass_1aa3fce0__1->GetZaxis()->SetLabelFont(42);
   frame_K_S0_mass_1aa3fce0__1->GetZaxis()->SetLabelSize(0.05);
   frame_K_S0_mass_1aa3fce0__1->GetZaxis()->SetTitleSize(0.05);
   frame_K_S0_mass_1aa3fce0__1->GetZaxis()->SetTitleOffset(1);
   frame_K_S0_mass_1aa3fce0__1->GetZaxis()->SetTitleFont(42);
   frame_K_S0_mass_1aa3fce0__1->Draw("FUNC");
   
   Double_t model_NormoBK_S0_masscB_CompoBGausscObackgroundcB_fx1[137] = { 0.397998, 0.398, 0.4, 0.402, 0.404, 0.406, 0.408, 0.41, 0.412, 0.414, 0.416, 0.418, 0.42, 0.422, 0.424, 0.426, 0.428,
   0.43, 0.432, 0.434, 0.436, 0.438, 0.44, 0.442, 0.444, 0.446, 0.448, 0.45, 0.452, 0.454, 0.456, 0.458, 0.46,
   0.462, 0.464, 0.466, 0.468, 0.469, 0.47, 0.471, 0.472, 0.473, 0.474, 0.475, 0.476, 0.477, 0.478, 0.479, 0.48,
   0.481, 0.482, 0.484, 0.485, 0.486, 0.4865, 0.487, 0.4875, 0.488, 0.4885, 0.489, 0.4895, 0.49, 0.4905, 0.491, 0.4915,
   0.492, 0.4925, 0.493, 0.4935, 0.494, 0.4945, 0.495, 0.4955, 0.496, 0.497, 0.498, 0.5, 0.501, 0.502, 0.503, 0.504,
   0.505, 0.506, 0.507, 0.508, 0.509, 0.51, 0.511, 0.512, 0.513, 0.514, 0.516, 0.518, 0.52, 0.522, 0.524, 0.526,
   0.528, 0.53, 0.532, 0.534, 0.536, 0.538, 0.54, 0.542, 0.544, 0.546, 0.548, 0.55, 0.552, 0.554, 0.556, 0.558,
   0.56, 0.562, 0.564, 0.566, 0.568, 0.57, 0.572, 0.574, 0.576, 0.578, 0.58, 0.582, 0.584, 0.586, 0.588, 0.59,
   0.592, 0.594, 0.596, 0.598, 0.6, 0.6, 0.602, 0.602002 };
   Double_t model_NormoBK_S0_masscB_CompoBGausscObackgroundcB_fy1[137] = { 0, 395.1721, 395.1721, 390.0397, 384.9377, 379.8662, 374.8254, 369.8153, 364.8363, 359.8885, 354.9719, 350.0869, 345.2335, 340.4119, 335.6222, 330.8647, 326.1394,
   321.4466, 316.7864, 312.159, 307.5645, 303.0032, 298.475, 293.9804, 289.5194, 285.0929, 280.7026, 276.3541, 272.0641, 267.8799, 263.9267, 260.515, 258.361,
   259.0004, 265.4821, 283.39, 322.1005, 353.5179, 395.926, 451.9419, 524.4417, 616.4559, 731.0149, 870.9445, 1038.615, 1235.656, 1462.657, 1718.873, 2001.975,
   2307.878, 2630.662, 3294.58, 3616.041, 3915.909, 4054.234, 4182.996, 4300.898, 4406.721, 4499.353, 4577.798, 4641.202, 4688.864, 4720.252, 4735.008, 4732.959,
   4714.117, 4678.681, 4627.033, 4559.729, 4477.493, 4381.198, 4271.856, 4150.597, 4018.65, 3727.98, 3410.857, 2741.732, 2410.335, 2092.711, 1795.555, 1523.746,
   1280.354, 1066.779, 882.9755, 727.7255, 598.9421, 493.9631, 409.8185, 343.4559, 291.9169, 252.4629, 200.3847, 171.7554, 156.2151, 147.4717, 142.0268, 138.0881,
   134.8061, 131.807, 128.9379, 126.1405, 123.394, 120.6918, 118.0317, 115.4133, 112.8367, 110.3018, 107.809, 105.3582, 102.9498, 100.5838, 98.26035, 95.97972,
   93.74203, 91.54743, 89.39609, 87.28816, 85.22382, 83.20321, 81.2265, 79.29386, 77.40544, 75.56141, 73.76192, 72.00714, 70.29723, 68.63236, 67.01267, 65.43834,
   63.90952, 62.42637, 60.98907, 59.59776, 58.25262, 58.25262, 58.25262, 0 };
   TGraph *graph = new TGraph(137,model_NormoBK_S0_masscB_CompoBGausscObackgroundcB_fx1,model_NormoBK_S0_masscB_CompoBGausscObackgroundcB_fy1);
   graph->SetName("model_Norm[K_S0_mass]_Comp[Gauss,background]");
   graph->SetTitle("Projection of model");

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#33ccff");
   graph->SetFillColor(ci);
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#33ccff");
   graph->SetLineColor(ci);
   graph->SetLineWidth(3);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.2);
   graph->Draw("f");
   
   Double_t h_K_S0_mass_fx3001[50] = { 0.402, 0.406, 0.41, 0.414, 0.418, 0.422, 0.426, 0.43, 0.434, 0.438, 0.442, 0.446, 0.45, 0.454, 0.458, 0.462, 0.466,
   0.47, 0.474, 0.478, 0.482, 0.486, 0.49, 0.494, 0.498, 0.502, 0.506, 0.51, 0.514, 0.518, 0.522, 0.526, 0.53,
   0.534, 0.538, 0.542, 0.546, 0.55, 0.554, 0.558, 0.562, 0.566, 0.57, 0.574, 0.578, 0.582, 0.586, 0.59, 0.594,
   0.598 };
   Double_t h_K_S0_mass_fy3001[50] = { 465, 407, 358, 371, 345, 346, 315, 293, 273, 255, 247, 256, 237, 245, 281, 293, 324,
   461, 689, 1218, 2507, 4308, 5150, 4281, 2998, 1933, 1119, 632, 395, 264, 199, 150, 143,
   105, 112, 92, 82, 81, 64, 79, 76, 78, 82, 72, 55, 68, 82, 80, 81,
   65 };
   Double_t h_K_S0_mass_felx3001[50] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0 };
   Double_t h_K_S0_mass_fely3001[50] = { 21.56386, 20.17424, 18.92089, 19.26136, 18.57418, 18.60108, 17.74824, 17.11724, 16.52271, 15.96872, 15.71623, 16, 15.3948, 15.65248, 16.76305, 17.11724, 18,
   21.47091, 26.24881, 34.89986, 50.06995, 65.63536, 71.7635, 65.42935, 54.754, 43.9659, 33.45146, 25.13961, 19.87461, 16.24808, 14.10674, 12.24745, 11.95826,
   10.24695, 10.58301, 9.591663, 9.055385, 9, 8, 8.888194, 8.717798, 8.831761, 9.055385, 8.485281, 7.416198, 8.246211, 9.055385, 8.944272, 9,
   8.062258 };
   Double_t h_K_S0_mass_fehx3001[50] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0 };
   Double_t h_K_S0_mass_fehy3001[50] = { 21.56386, 20.17424, 18.92089, 19.26136, 18.57418, 18.60108, 17.74824, 17.11724, 16.52271, 15.96872, 15.71623, 16, 15.3948, 15.65248, 16.76305, 17.11724, 18,
   21.47091, 26.24881, 34.89986, 50.06995, 65.63536, 71.7635, 65.42935, 54.754, 43.9659, 33.45146, 25.13961, 19.87461, 16.24808, 14.10674, 12.24745, 11.95826,
   10.24695, 10.58301, 9.591663, 9.055385, 9, 8, 8.888194, 8.717798, 8.831761, 9.055385, 8.485281, 7.416198, 8.246211, 9.055385, 8.944272, 9,
   8.062258 };
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(50,h_K_S0_mass_fx3001,h_K_S0_mass_fy3001,h_K_S0_mass_felx3001,h_K_S0_mass_fehx3001,h_K_S0_mass_fely3001,h_K_S0_mass_fehy3001);
   grae->SetName("h_K_S0_mass");
   grae->SetTitle("Histogram of K_S0_mass_plot__K_S0_mass");
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(8);
   grae->SetMarkerSize(1.2);
   grae->Draw("p");
   
   Double_t h_K_S0_mass_fx3002[50] = { 0.402, 0.406, 0.41, 0.414, 0.418, 0.422, 0.426, 0.43, 0.434, 0.438, 0.442, 0.446, 0.45, 0.454, 0.458, 0.462, 0.466,
   0.47, 0.474, 0.478, 0.482, 0.486, 0.49, 0.494, 0.498, 0.502, 0.506, 0.51, 0.514, 0.518, 0.522, 0.526, 0.53,
   0.534, 0.538, 0.542, 0.546, 0.55, 0.554, 0.558, 0.562, 0.566, 0.57, 0.574, 0.578, 0.582, 0.586, 0.59, 0.594,
   0.598 };
   Double_t h_K_S0_mass_fy3002[50] = { 465, 407, 358, 371, 345, 346, 315, 293, 273, 255, 247, 256, 237, 245, 281, 293, 324,
   461, 689, 1218, 2507, 4308, 5150, 4281, 2998, 1933, 1119, 632, 395, 264, 199, 150, 143,
   105, 112, 92, 82, 81, 64, 79, 76, 78, 82, 72, 55, 68, 82, 80, 81,
   65 };
   Double_t h_K_S0_mass_felx3002[50] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0 };
   Double_t h_K_S0_mass_fely3002[50] = { 21.56386, 20.17424, 18.92089, 19.26136, 18.57418, 18.60108, 17.74824, 17.11724, 16.52271, 15.96872, 15.71623, 16, 15.3948, 15.65248, 16.76305, 17.11724, 18,
   21.47091, 26.24881, 34.89986, 50.06995, 65.63536, 71.7635, 65.42935, 54.754, 43.9659, 33.45146, 25.13961, 19.87461, 16.24808, 14.10674, 12.24745, 11.95826,
   10.24695, 10.58301, 9.591663, 9.055385, 9, 8, 8.888194, 8.717798, 8.831761, 9.055385, 8.485281, 7.416198, 8.246211, 9.055385, 8.944272, 9,
   8.062258 };
   Double_t h_K_S0_mass_fehx3002[50] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0 };
   Double_t h_K_S0_mass_fehy3002[50] = { 21.56386, 20.17424, 18.92089, 19.26136, 18.57418, 18.60108, 17.74824, 17.11724, 16.52271, 15.96872, 15.71623, 16, 15.3948, 15.65248, 16.76305, 17.11724, 18,
   21.47091, 26.24881, 34.89986, 50.06995, 65.63536, 71.7635, 65.42935, 54.754, 43.9659, 33.45146, 25.13961, 19.87461, 16.24808, 14.10674, 12.24745, 11.95826,
   10.24695, 10.58301, 9.591663, 9.055385, 9, 8, 8.888194, 8.717798, 8.831761, 9.055385, 8.485281, 7.416198, 8.246211, 9.055385, 8.944272, 9,
   8.062258 };
   grae = new TGraphAsymmErrors(50,h_K_S0_mass_fx3002,h_K_S0_mass_fy3002,h_K_S0_mass_felx3002,h_K_S0_mass_fehx3002,h_K_S0_mass_fely3002,h_K_S0_mass_fehy3002);
   grae->SetName("h_K_S0_mass");
   grae->SetTitle("Histogram of K_S0_mass_plot__K_S0_mass");
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(8);
   grae->SetMarkerSize(1.2);
   grae->Draw("p");
   
   Double_t model_NormoBK_S0_masscB_CompoBbackgroundcB_fx2[106] = { 0.397998, 0.398, 0.4, 0.402, 0.404, 0.406, 0.408, 0.41, 0.412, 0.414, 0.416, 0.418, 0.42, 0.422, 0.424, 0.426, 0.428,
   0.43, 0.432, 0.434, 0.436, 0.438, 0.44, 0.442, 0.444, 0.446, 0.448, 0.45, 0.452, 0.454, 0.456, 0.458, 0.46,
   0.462, 0.464, 0.466, 0.468, 0.47, 0.472, 0.474, 0.476, 0.478, 0.48, 0.482, 0.484, 0.486, 0.488, 0.49, 0.492,
   0.494, 0.496, 0.498, 0.5, 0.502, 0.504, 0.506, 0.508, 0.51, 0.512, 0.514, 0.516, 0.518, 0.52, 0.522, 0.524,
   0.526, 0.528, 0.53, 0.532, 0.534, 0.536, 0.538, 0.54, 0.542, 0.544, 0.546, 0.548, 0.55, 0.552, 0.554, 0.556,
   0.558, 0.56, 0.562, 0.564, 0.566, 0.568, 0.57, 0.572, 0.574, 0.576, 0.578, 0.58, 0.582, 0.584, 0.586, 0.588,
   0.59, 0.592, 0.594, 0.596, 0.598, 0.6, 0.6, 0.602, 0.602002 };
   Double_t model_NormoBK_S0_masscB_CompoBbackgroundcB_fy2[106] = { 0, 395.1721, 395.1721, 390.0397, 384.9377, 379.8662, 374.8254, 369.8153, 364.8363, 359.8885, 354.9719, 350.0869, 345.2335, 340.4119, 335.6222, 330.8647, 326.1394,
   321.4466, 316.7864, 312.159, 307.5645, 303.0031, 298.475, 293.9803, 289.5191, 285.0917, 280.6982, 276.3387, 272.0135, 267.7227, 263.4664, 259.2448, 255.058,
   250.9063, 246.7898, 242.7086, 238.663, 234.653, 230.6789, 226.7407, 222.8387, 218.9731, 215.1439, 211.3513, 207.5956, 203.8768, 200.1951, 196.5508, 192.9438,
   189.3745, 185.843, 182.3494, 178.8939, 175.4766, 172.0978, 168.7575, 165.456, 162.1933, 158.9698, 155.7854, 152.6405, 149.535, 146.4693, 143.4435, 140.4577,
   137.5121, 134.6068, 131.742, 128.918, 126.1347, 123.3925, 120.6914, 118.0316, 115.4133, 112.8367, 110.3018, 107.809, 105.3582, 102.9498, 100.5838, 98.26035,
   95.97972, 93.74203, 91.54743, 89.39609, 87.28816, 85.22382, 83.20321, 81.2265, 79.29386, 77.40544, 75.56141, 73.76192, 72.00714, 70.29723, 68.63236, 67.01267,
   65.43834, 63.90952, 62.42637, 60.98907, 59.59776, 58.25262, 58.25262, 58.25262, 0 };
   graph = new TGraph(106,model_NormoBK_S0_masscB_CompoBbackgroundcB_fx2,model_NormoBK_S0_masscB_CompoBbackgroundcB_fy2);
   graph->SetName("model_Norm[K_S0_mass]_Comp[background]");
   graph->SetTitle("Projection of model");

   ci = TColor::GetColor("#cccccc");
   graph->SetFillColor(ci);
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#cccccc");
   graph->SetLineColor(ci);
   graph->SetLineWidth(3);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.2);
   graph->Draw("f");
   
   Double_t model_NormoBK_S0_masscB_fx3[137] = { 0.397998, 0.398, 0.4, 0.402, 0.404, 0.406, 0.408, 0.41, 0.412, 0.414, 0.416, 0.418, 0.42, 0.422, 0.424, 0.426, 0.428,
   0.43, 0.432, 0.434, 0.436, 0.438, 0.44, 0.442, 0.444, 0.446, 0.448, 0.45, 0.452, 0.454, 0.456, 0.458, 0.46,
   0.462, 0.464, 0.466, 0.468, 0.469, 0.47, 0.471, 0.472, 0.473, 0.474, 0.475, 0.476, 0.477, 0.478, 0.479, 0.48,
   0.481, 0.482, 0.484, 0.485, 0.486, 0.4865, 0.487, 0.4875, 0.488, 0.4885, 0.489, 0.4895, 0.49, 0.4905, 0.491, 0.4915,
   0.492, 0.4925, 0.493, 0.4935, 0.494, 0.4945, 0.495, 0.4955, 0.496, 0.497, 0.498, 0.5, 0.501, 0.502, 0.503, 0.504,
   0.505, 0.506, 0.507, 0.508, 0.509, 0.51, 0.511, 0.512, 0.513, 0.514, 0.516, 0.518, 0.52, 0.522, 0.524, 0.526,
   0.528, 0.53, 0.532, 0.534, 0.536, 0.538, 0.54, 0.542, 0.544, 0.546, 0.548, 0.55, 0.552, 0.554, 0.556, 0.558,
   0.56, 0.562, 0.564, 0.566, 0.568, 0.57, 0.572, 0.574, 0.576, 0.578, 0.58, 0.582, 0.584, 0.586, 0.588, 0.59,
   0.592, 0.594, 0.596, 0.598, 0.6, 0.6, 0.602, 0.602002 };
   Double_t model_NormoBK_S0_masscB_fy3[137] = { 0, 395.1721, 395.1721, 390.0397, 384.9377, 379.8662, 374.8254, 369.8153, 364.8363, 359.8885, 354.9719, 350.0869, 345.2335, 340.4119, 335.6222, 330.8647, 326.1394,
   321.4466, 316.7864, 312.159, 307.5645, 303.0032, 298.475, 293.9804, 289.5194, 285.0929, 280.7026, 276.3541, 272.0641, 267.8799, 263.9267, 260.515, 258.361,
   259.0004, 265.4821, 283.39, 322.1005, 353.5179, 395.926, 451.9419, 524.4417, 616.4559, 731.0149, 870.9445, 1038.615, 1235.656, 1462.657, 1718.873, 2001.975,
   2307.878, 2630.662, 3294.58, 3616.041, 3915.909, 4054.234, 4182.996, 4300.898, 4406.721, 4499.353, 4577.798, 4641.202, 4688.864, 4720.252, 4735.008, 4732.959,
   4714.117, 4678.681, 4627.033, 4559.729, 4477.493, 4381.198, 4271.856, 4150.597, 4018.65, 3727.98, 3410.857, 2741.732, 2410.335, 2092.711, 1795.555, 1523.746,
   1280.354, 1066.779, 882.9755, 727.7255, 598.9421, 493.9631, 409.8185, 343.4559, 291.9169, 252.4629, 200.3847, 171.7554, 156.2151, 147.4717, 142.0268, 138.0881,
   134.8061, 131.807, 128.9379, 126.1405, 123.394, 120.6918, 118.0317, 115.4133, 112.8367, 110.3018, 107.809, 105.3582, 102.9498, 100.5838, 98.26035, 95.97972,
   93.74203, 91.54743, 89.39609, 87.28816, 85.22382, 83.20321, 81.2265, 79.29386, 77.40544, 75.56141, 73.76192, 72.00714, 70.29723, 68.63236, 67.01267, 65.43834,
   63.90952, 62.42637, 60.98907, 59.59776, 58.25262, 58.25262, 58.25262, 0 };
   graph = new TGraph(137,model_NormoBK_S0_masscB_fx3,model_NormoBK_S0_masscB_fy3);
   graph->SetName("model_Norm[K_S0_mass]");
   graph->SetTitle("Projection of model");
   graph->SetFillStyle(1000);
   graph->SetLineWidth(3);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.2);
   graph->Draw("l");
   
   Double_t h_K_S0_mass_fx3003[50] = { 0.402, 0.406, 0.41, 0.414, 0.418, 0.422, 0.426, 0.43, 0.434, 0.438, 0.442, 0.446, 0.45, 0.454, 0.458, 0.462, 0.466,
   0.47, 0.474, 0.478, 0.482, 0.486, 0.49, 0.494, 0.498, 0.502, 0.506, 0.51, 0.514, 0.518, 0.522, 0.526, 0.53,
   0.534, 0.538, 0.542, 0.546, 0.55, 0.554, 0.558, 0.562, 0.566, 0.57, 0.574, 0.578, 0.582, 0.586, 0.59, 0.594,
   0.598 };
   Double_t h_K_S0_mass_fy3003[50] = { 465, 407, 358, 371, 345, 346, 315, 293, 273, 255, 247, 256, 237, 245, 281, 293, 324,
   461, 689, 1218, 2507, 4308, 5150, 4281, 2998, 1933, 1119, 632, 395, 264, 199, 150, 143,
   105, 112, 92, 82, 81, 64, 79, 76, 78, 82, 72, 55, 68, 82, 80, 81,
   65 };
   Double_t h_K_S0_mass_felx3003[50] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0 };
   Double_t h_K_S0_mass_fely3003[50] = { 21.56386, 20.17424, 18.92089, 19.26136, 18.57418, 18.60108, 17.74824, 17.11724, 16.52271, 15.96872, 15.71623, 16, 15.3948, 15.65248, 16.76305, 17.11724, 18,
   21.47091, 26.24881, 34.89986, 50.06995, 65.63536, 71.7635, 65.42935, 54.754, 43.9659, 33.45146, 25.13961, 19.87461, 16.24808, 14.10674, 12.24745, 11.95826,
   10.24695, 10.58301, 9.591663, 9.055385, 9, 8, 8.888194, 8.717798, 8.831761, 9.055385, 8.485281, 7.416198, 8.246211, 9.055385, 8.944272, 9,
   8.062258 };
   Double_t h_K_S0_mass_fehx3003[50] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0 };
   Double_t h_K_S0_mass_fehy3003[50] = { 21.56386, 20.17424, 18.92089, 19.26136, 18.57418, 18.60108, 17.74824, 17.11724, 16.52271, 15.96872, 15.71623, 16, 15.3948, 15.65248, 16.76305, 17.11724, 18,
   21.47091, 26.24881, 34.89986, 50.06995, 65.63536, 71.7635, 65.42935, 54.754, 43.9659, 33.45146, 25.13961, 19.87461, 16.24808, 14.10674, 12.24745, 11.95826,
   10.24695, 10.58301, 9.591663, 9.055385, 9, 8, 8.888194, 8.717798, 8.831761, 9.055385, 8.485281, 7.416198, 8.246211, 9.055385, 8.944272, 9,
   8.062258 };
   grae = new TGraphAsymmErrors(50,h_K_S0_mass_fx3003,h_K_S0_mass_fy3003,h_K_S0_mass_felx3003,h_K_S0_mass_fehx3003,h_K_S0_mass_fely3003,h_K_S0_mass_fehy3003);
   grae->SetName("h_K_S0_mass");
   grae->SetTitle("Histogram of K_S0_mass_plot__K_S0_mass");
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(8);
   grae->SetMarkerSize(1.2);
   grae->Draw("pe1");
   
   TH1D *frame_K_S0_mass_1aa3fce0__2 = new TH1D("frame_K_S0_mass_1aa3fce0__2","A RooPlot of \"mass\"",100,0.4,0.6);
   frame_K_S0_mass_1aa3fce0__2->SetBinContent(1,5482.852);
   frame_K_S0_mass_1aa3fce0__2->SetMaximum(5482.852);
   frame_K_S0_mass_1aa3fce0__2->SetEntries(1);
   frame_K_S0_mass_1aa3fce0__2->SetDirectory(nullptr);
   frame_K_S0_mass_1aa3fce0__2->SetStats(0);
   frame_K_S0_mass_1aa3fce0__2->SetMarkerStyle(4);
   frame_K_S0_mass_1aa3fce0__2->SetMarkerSize(0.02);
   frame_K_S0_mass_1aa3fce0__2->GetXaxis()->SetTitle("m(#pi^{+}#pi^{-}) [GeV]");
   frame_K_S0_mass_1aa3fce0__2->GetXaxis()->SetNdivisions(605);
   frame_K_S0_mass_1aa3fce0__2->GetXaxis()->SetLabelFont(42);
   frame_K_S0_mass_1aa3fce0__2->GetXaxis()->SetLabelSize(0.05);
   frame_K_S0_mass_1aa3fce0__2->GetXaxis()->SetTitleSize(0.05);
   frame_K_S0_mass_1aa3fce0__2->GetXaxis()->SetTitleOffset(1.4);
   frame_K_S0_mass_1aa3fce0__2->GetXaxis()->SetTitleFont(42);
   frame_K_S0_mass_1aa3fce0__2->GetYaxis()->SetTitle("Candidates / (4 MeV)");
   frame_K_S0_mass_1aa3fce0__2->GetYaxis()->SetLabelFont(42);
   frame_K_S0_mass_1aa3fce0__2->GetYaxis()->SetLabelSize(0.05);
   frame_K_S0_mass_1aa3fce0__2->GetYaxis()->SetTitleSize(0.05);
   frame_K_S0_mass_1aa3fce0__2->GetYaxis()->SetTitleOffset(1.4);
   frame_K_S0_mass_1aa3fce0__2->GetYaxis()->SetTitleFont(42);
   frame_K_S0_mass_1aa3fce0__2->GetZaxis()->SetLabelFont(42);
   frame_K_S0_mass_1aa3fce0__2->GetZaxis()->SetLabelSize(0.05);
   frame_K_S0_mass_1aa3fce0__2->GetZaxis()->SetTitleSize(0.05);
   frame_K_S0_mass_1aa3fce0__2->GetZaxis()->SetTitleOffset(1);
   frame_K_S0_mass_1aa3fce0__2->GetZaxis()->SetTitleFont(42);
   frame_K_S0_mass_1aa3fce0__2->Draw("AXISSAME");
   
   TPaveText *pt = new TPaveText(0,0,0,0,"brNDC");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   TText *pt_LaTex = pt->AddText("#it{#bf{sPHENIX}} Internal");
   pt_LaTex = pt->AddText("#it{p}+#it{p} #sqrt{s} = 200 GeV");
   pt_LaTex = pt->AddText("#mu = 491#kern[-0.4]{ #pm} 0 MeV");
   pt_LaTex = pt->AddText("#sigma = 8.2#kern[-0.4]{ #pm} 0.1 MeV");
   pt_LaTex = pt->AddText("Yield = 23371#kern[-0.3]{ #pm} 166");
   pt->Draw();
   
   pt = new TPaveText(0,0,0,0,"brNDC");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->SetTextSize(0.05);
   pt_LaTex = pt->AddText("05/22/2025");
   pt->Draw();
   
   TLegend *leg = new TLegend(0,0,0,0,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("h_K_S0_mass","Data","PE2");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("model_Norm[K_S0_mass]","Fit","L");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("model_Norm[K_S0_mass]_Comp[Gauss,background]","K_{S}^{0}#rightarrow#pi^{+}#pi^{-}","f");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("model_Norm[K_S0_mass]_Comp[background]","Comb. Bkg.","f");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   massFitCanvas->Modified();
   massFitCanvas->SetSelected(massFitCanvas);
}
