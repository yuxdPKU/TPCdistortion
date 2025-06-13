#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>

void plot_clusters() {
    // 打开ROOT文件
    TFile* file = TFile::Open("root/Reconstructed/53877/clusters_seeds_53877-0.root_resid.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file!" << std::endl;
        return;
    }

    // 获取residualtree
    TTree* tree = (TTree*)file->Get("residualtree");
    if (!tree) {
        std::cerr << "Error: Could not find residualtree!" << std::endl;
        file->Close();
        return;
    }

    // 创建画布
    TCanvas* canvas = new TCanvas("canvas", "Cluster Positions", 800, 600);
    
    // 设置样式
    gStyle->SetMarkerSize(0.2);
    gStyle->SetMarkerStyle(20);
    
    // 创建TGraph
    int n = tree->Draw("clusgy:clusgx", "cluslayer>=3 && cluslayer<=6", "goff");
    TGraph* graph = new TGraph(n, tree->GetV2(), tree->GetV1());
    
    // 设置图形属性
    graph->SetTitle("Cluster Positions;clusgx;clusgy");
    graph->Draw("AP");
    
    // 保存图像
    canvas->SaveAs("cluster_positions.pdf");
    
    // 清理
    delete graph;
    delete canvas;
    file->Close();
}
