#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BUFSIZE 2000
#define MAXSEQ 1000
#define LOG_BASE 2.0

char *labels[4] = {"A", "C", "G", "T"};

typedef struct {
    char name[BUFSIZE];
    char seq[BUFSIZE];
} SEQ;

SEQ g_motif[MAXSEQ];
SEQ g_pro[MAXSEQ];

// A:0, C:1, G:2, T:3
int base2num(char base){
    switch(base){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

// モチーフ配列読み込み
int read_multi_seq(char* filename){
    FILE *fp = fopen(filename, "r");
    if(fp == NULL){
        printf("motif_file open error.\n");
        exit(1);
    }

    char line[BUFSIZE];
    int count = 0;
    while(fgets(line, BUFSIZE, fp) != NULL){
        if(line[0] == '>'){
            strcpy(g_motif[count].name, line + 1);
            g_motif[count].name[strcspn(g_motif[count].name, "\n")] = '\0';
        } else {
            line[strcspn(line, "\n")] = '\0';
            strcpy(g_motif[count].seq, line);
            count++;
        }
    }
    fclose(fp);
    return count;
}

// プロモータ配列読み込み
int read_promoter(char* filename){
    FILE *fp = fopen(filename, "r");
    if(fp == NULL){
        printf("promoter_file open error.\n");
        exit(1);
    }

    char line[BUFSIZE];
    int count = 0;
    while(fgets(line, BUFSIZE, fp) != NULL){
        if(line[0] == '>'){
            strcpy(g_pro[count].name, line + 1);
            g_pro[count].name[strcspn(g_pro[count].name, "\n")] = '\0';
        } else {
            line[strcspn(line, "\n")] = '\0';
            strcpy(g_pro[count].seq, line);
            count++;
        }
    }
    fclose(fp);
    return count;
}

// 頻度行列を作成
void make_freq_table(int freqtable[4][BUFSIZE], int len, int seq_num){
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < len; j++)
            freqtable[i][j] = 0;

    for(int i = 0; i < seq_num; i++){
        for(int j = 0; j < len; j++){
            int idx = base2num(g_motif[i].seq[j]);
            if(idx >= 0) freqtable[idx][j]++;
        }
    }
}

// log-oddsスコアを計算
void make_odds_score(float q[4][BUFSIZE], int freqtable[4][BUFSIZE], int len, int seq_num){
    float bg[4];
    float total_bg = 7519429 + 4637676 + 4637676 + 7519429;
    bg[0] = 7519429.0 / total_bg; // A
    bg[1] = 4637676.0 / total_bg; // C
    bg[2] = 4637676.0 / total_bg; // G
    bg[3] = 7519429.0 / total_bg; // T

    for(int j = 0; j < len; j++){
        for(int i = 0; i < 4; i++){
            float pij = (freqtable[i][j] + 1.0) / (seq_num + 4.0); // 疑似頻度1の加算平滑化
            q[i][j] = log(pij / bg[i]);
        }
    }
}

// スコア出力（スコア5以上のすべて，なければ最大スコア1件）
// 各プロモータごとに：スコア5.0以上をすべて表示，なければ最大スコア1件のみ表示
void score_promoters(int gene_num, int motif_len, float q[4][BUFSIZE]){
    printf("\n--- プロモータスコア（各プロモータごとに：スコア5.0以上すべて，なければ最大1件）---\n");

    for(int i = 0; i < gene_num; i++){
        int seq_len = strlen(g_pro[i].seq);
        float max_score = -INFINITY;
        int max_pos = -1;

        // スコア5以上のものを一時保存
        typedef struct {
            int pos;
            float score;
        } ScoreEntry;

        ScoreEntry entries[BUFSIZE];
        int entry_count = 0;

        for(int j = 0; j <= seq_len - motif_len; j++){
            float score = 0.0;
            for(int k = 0; k < motif_len; k++){
                int idx = base2num(g_pro[i].seq[j + k]);
                if(idx < 0){
                    score = -INFINITY;
                    break;
                }
                score += q[idx][k];
            }

            if(score >= 5.0){
                entries[entry_count].pos = j;
                entries[entry_count].score = score;
                entry_count++;
            }

            if(score > max_score){
                max_score = score;
                max_pos = j;
            }
        }

        if(entry_count > 0){
            for(int j = 0; j < entry_count; j++){
                printf("%s\tpos: %d\tscore: %.2f\n", g_pro[i].name, entries[j].pos, entries[j].score);
            }
        } else {
            printf("%s\tpos: %d\tscore: %.2f\n", g_pro[i].name, max_pos, max_score);
        }
    }
}




int main(int argc, char* argv[]){
    if(argc < 3){
        printf("Usage: %s motif_file promoter_file\n", argv[0]);
        return 1;
    }

    int seq_num = read_multi_seq(argv[1]);
    int len = strlen(g_motif[0].seq);

    int freqtable[4][BUFSIZE];
    float q[4][BUFSIZE];

    make_freq_table(freqtable, len, seq_num);
    make_odds_score(q, freqtable, len, seq_num);

    // 頻度出力
    printf("列ごとの塩基頻度：\n");
    for(int i = 0; i < 4; i++){
        printf("%s: ", labels[i]);
        for(int j = 0; j < len; j++){
            printf("%2d ", freqtable[i][j]);
        }
        printf("\n");
    }

    // ログオッズスコア出力
    printf("\n列ごとのオッズスコア（log-odds）：\n");
    for(int i = 0; i < 4; i++){
        printf("%s: ", labels[i]);
        for(int j = 0; j < len; j++){
            printf("%5.2f ", q[i][j]);
        }
        printf("\n");
    }

    int gene_num = read_promoter(argv[2]);
    score_promoters(gene_num, len, q);

    return 0;
}
