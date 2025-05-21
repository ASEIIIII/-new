#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXLEN 2000
#define MAXSEQ 1000

char *base_labels[4] = {"A", "C", "G", "T"};

typedef struct {
    char name[MAXLEN];
    char seq[MAXLEN];
} Sequence;

Sequence motifs[MAXSEQ];
Sequence promoters[MAXSEQ];

int base_to_index(char base){
    switch(base){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

int load_motifs(char* filename){
    FILE *fp = fopen(filename, "r");
    if(fp == NULL){
        perror("モチーフファイル開けなかった");
        exit(1);
    }

    char line[MAXLEN];
    int count = 0;
    while(fgets(line, MAXLEN, fp) != NULL){
        if(line[0] == '>'){
            strcpy(motifs[count].name, line + 1);
            motifs[count].name[strcspn(motifs[count].name, "\n")] = '\0';
        } else {
            line[strcspn(line, "\n")] = '\0';
            strcpy(motifs[count].seq, line);
            count++;
        }
    }
    fclose(fp);
    return count;
}

int load_promoters(char* filename){
    FILE *fp = fopen(filename, "r");
    if(fp == NULL){
        perror("プロモータファイル開けなかった");
        exit(1);
    }

    char line[MAXLEN];
    int count = 0;
    while(fgets(line, MAXLEN, fp) != NULL){
        if(line[0] == '>'){
            strcpy(promoters[count].name, line + 1);
            promoters[count].name[strcspn(promoters[count].name, "\n")] = '\0';
        } else {
            line[strcspn(line, "\n")] = '\0';
            strcpy(promoters[count].seq, line);
            count++;
        }
    }
    fclose(fp);
    return count;
}

void calc_freq_table(int freq[4][MAXLEN], int motif_len, int motif_num){
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < motif_len; j++)
            freq[i][j] = 0;

    for(int i = 0; i < motif_num; i++){
        for(int j = 0; j < motif_len; j++){
            int idx = base_to_index(motifs[i].seq[j]);
            if(idx >= 0) freq[idx][j]++;
        }
    }
}

void calc_score_matrix(float score[4][MAXLEN], int freq[4][MAXLEN], int motif_len, int motif_num){
    float bg[4];
    float total = 7519429 + 4637676 + 4637676 + 7519429;
    bg[0] = 7519429.0 / total; 
    bg[1] = 4637676.0 / total; 
    bg[2] = 4637676.0 / total; 
    bg[3] = 7519429.0 / total; 

    for(int j = 0; j < motif_len; j++){
        for(int i = 0; i < 4; i++){
            float prob = (freq[i][j] + 1.0) / (motif_num + 4.0); 
            score[i][j] = log(prob / bg[i]);
        }
    }
}

void search_binding_sites(int promoter_count, int motif_len, float score_matrix[4][MAXLEN]){
    for(int i = 0; i < promoter_count; i++){
        char* pseq = promoters[i].seq;
        int plen = strlen(pseq);
        float best_score = -INFINITY;
        int best_pos = -1;

        printf("\n%s\n", promoters[i].name);

        for(int j = 0; j <= plen - motif_len; j++){
            float score = 0;
            for(int k = 0; k < motif_len; k++){
                int b = base_to_index(pseq[j + k]);
                if(b < 0){
                    score = -INFINITY;
                    break;
                }
                score += score_matrix[b][k];
            }

            if(score >= 5.0){
                printf("  * ヒット位置: %d, スコア: %.2f\n", j, score);
            }

            if(score > best_score){
                best_score = score;
                best_pos = j;
            }
        }

        if(best_score < 5.0){
            printf("  最も高いスコア: %.2f (位置: %d)\n", best_score, best_pos);
        }
    }
}

int main(int argc, char* argv[]){

    int motif_num = load_motifs(argv[1]);
    int motif_len = strlen(motifs[0].seq);
    
    int freq[4][MAXLEN];
    float score[4][MAXLEN];

    calc_freq_table(freq, motif_len, motif_num);
    calc_score_matrix(score, freq, motif_len, motif_num);

    int promoter_count = load_promoters(argv[2]);
    search_binding_sites(promoter_count, motif_len, score);

    return 0;
}
