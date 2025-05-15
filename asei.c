#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BUFSIZE 1024
#define MAX_SEQ_NUM 30
#define MAX_GENE_NUM 8

char g_motif[MAX_SEQ_NUM][BUFSIZE];

struct promoter {
    char name[BUFSIZE];
    char seq[BUFSIZE];
} g_pro[MAX_GENE_NUM];

int read_multi_seq(char* filename){
    int seq_num = 0;
    char buffer[BUFSIZE];
    FILE *fp = fopen(filename,"r");

    if(fp == NULL){
        printf("motif_region_file open error.\n");
        exit(1);
    }

    while(fscanf(fp, "%s", buffer) != EOF){
        strcpy(g_motif[seq_num], buffer);
        seq_num++;
    }
    fclose(fp);
    return seq_num;
}

int main(int argc, char* argv[]){
    if(argc < 2){
        printf("Usage: %s motif_file\n", argv[0]);
        return 1;
    }

    int seq_num = read_multi_seq(argv[1]);
    int len = strlen(g_motif[0]);

    int freqtable[4][BUFSIZE] = {0};  // A, C, G, T

    // 塩基頻度カウント
    for(int i = 0; i < seq_num; i++){
        for(int j = 0; j < len; j++){
            switch(g_motif[i][j]){
                case 'A': freqtable[0][j]++; break;
                case 'C': freqtable[1][j]++; break;
                case 'G': freqtable[2][j]++; break;
                case 'T': freqtable[3][j]++; break;
            }
        }
    }

    // 頻度表出力
    const char *labels[4] = {"A", "C", "G", "T"};
    printf("列ごとの塩基頻度：\n");
    for(int i = 0; i < 4; i++){
        printf("%s: ", labels[i]);
        for(int j = 0; j < len; j++){
            printf("%2d ", freqtable[i][j]);
        }
        printf("\n");
    }

    // 背景確率
    float bg_A = 7519429.0 / (7519429 + 4637676 + 4637676 + 7519429);
    float bg_C = 4637676.0 / (7519429 + 4637676 + 4637676 + 7519429);
    float bg_G = 4637676.0 / (7519429 + 4637676 + 4637676 + 7519429);
    float bg_T = 7519429.0 / (7519429 + 4637676 + 4637676 + 7519429);

    float bg[4] = {bg_A, bg_C, bg_G, bg_T};
    float q[4][BUFSIZE] = {0};

    // オッズスコア計算（加算平滑化付き）
    for(int j = 0; j < len; j++){
        float total = freqtable[0][j] + freqtable[1][j] + freqtable[2][j] + freqtable[3][j] + 4; // 加算平滑化
        for(int i = 0; i < 4; i++){
            q[i][j] = log((freqtable[i][j] + 1.0) / total / bg[i]); // log-odds
        }
    }

    // オッズスコア出力
    printf("\n列ごとのオッズスコア（log-odds）:\n");
    for(int i = 0; i < 4; i++){
        printf("%s: ", labels[i]);
        for(int j = 0; j < len; j++){
            printf("%5.2f ", q[i][j]);
        }
        printf("\n");
    }
    return 0;
}

    
