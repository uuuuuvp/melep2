#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "melpe.h"
#include "npp.h"
#include "sc1200.h"

// 1200
#define FRAME_SIZE 180 // 将每帧的样本数更改为 180
#define BITS_PER_FRAME 54 // 每帧的比特数，根据帧长计算
#define SAMPLE_RATE 8000 // 采样率为 8kHz
#define BITS_PER_SAMPLE 16

// struct melp_param pa;
// struct melp_param pat;

typedef struct {
    char RIFF[4];// RIFF标识，表示这是一个RIFF文件
    uint32_t fileSize;// 文件总大小，包括文件头和音频数据
    char WAVE[4];// WAVE标识，表示文件格式为WAVE
    char fmt[4];// "fmt "子块标识，描述音频格式的子块
    uint32_t fmtSize;
    uint16_t format;// 音频数据的格式
    uint16_t channels;// 音频的通道
    uint32_t sampleRate;// 音频的采样率
    uint32_t byteRate;// 字节速率
    uint16_t blockAlign;
    uint16_t bitsPerSample;
    char data[4];
    uint32_t dataSize;// 音频数据的大小
} WAVHeader;

// 生成 WAV 文件头的函数
void createWAVHeader(WAVHeader *header, uint32_t dataSize) {
    memcpy(header->RIFF, "RIFF", 4);
    header->fileSize = 36 + dataSize;
    memcpy(header->WAVE, "WAVE", 4);
    memcpy(header->fmt, "fmt ", 4);
    header->fmtSize = 16;
    header->format = 1; // PCM
    header->channels = 1;
    header->sampleRate = SAMPLE_RATE;
    header->byteRate = SAMPLE_RATE * header->channels * (BITS_PER_SAMPLE / 8);
    header->blockAlign = header->channels * (BITS_PER_SAMPLE / 8);
    header->bitsPerSample = BITS_PER_SAMPLE;
    memcpy(header->data, "data", 4);
    header->dataSize = dataSize;
}

// 主函数
int main() {
    FILE *input_file, *output_file; // 声明两个文件指针，用于输入和输出文件
    short speech[FRAME_SIZE]; // 存储每帧音频样本
    unsigned char bitstream[(BITS_PER_FRAME + 2) / 8]; // 存储每帧的比特流
    // unsigned char bitstream[(BITS_PER_FRAME + 2) / 8]; // 存储每帧的比特流2400
    WAVHeader header; // 定义wav文件头结构体实例

    input_file = fopen("input8KHz.wav", "rb");
    if (input_file == NULL) {
        printf("Error opening input file.\n");
        return 1;
    }
    output_file = fopen("action_14.wav", "wb");
    if (output_file == NULL) {
        printf("Error opening output file.\n");
        fclose(input_file);
        return 1;
    }

    melpe_i(); // 初始化编解码器

    // 使用fseek定位到输出文件，为了先占位写入 WAV 头的大小
    fseek(output_file, sizeof(WAVHeader), SEEK_SET); 

    uint32_t dataSize = 0; // 初始化音频数据部分的总大小为 0
    int i = 0;
    while (fread(speech, sizeof(short), FRAME_SIZE, input_file) == FRAME_SIZE) {
        // melpe_a(bitstream, speech); 
        // melpe_s(speech, bitstream); 
        // npp(speech, speech);
        analysis(speech, &pa);
        printf("%d",i+1);
             printf("pitch:%d\n",pa.pitch);
             printf("jitter:%d\n",pa.jitter);
             printf("uv_flag:%d",pa.uv_flag);
             printf("\nlsf:");
             for (int j = 0; j < LPC_ORD ; ++j) {
                  printf("%d ", pa.lsf[j]);
             }
             printf("\ngain:");
             for (int k = 0; k < NUM_GAINFR ; ++k) {
                  printf("%d ", pa.gain[k]);
             }
             printf("\nbpvc:");
             for (int l = 0; l < NUM_BANDS ; ++l) {
                  printf("%d ", pa.bpvc[l]);
             }
             printf("\nfs_mag:");
             for (int m = 0; m < NUM_HARM ; ++m) {
                  printf("%d ", pa.fs_mag[m]);
             }
             printf("\n\n");
             i++;
        // rmquan(&pat, speech, &pa);
        fwrite(speech, sizeof(short), FRAME_SIZE, output_file);  // 音频数据写入
        dataSize += FRAME_SIZE * sizeof(short); // 更新音频数据大小

    }
    // 使用 fseek 定位到输出文件的开头
    fseek(output_file, 0, SEEK_SET);
    // 并调用 createWAVHeader 函数生成 WAV 头
    createWAVHeader(&header, dataSize);
    // 并将其写入输出文件
    fwrite(&header, sizeof(WAVHeader), 1, output_file);

    fclose(input_file);
    fclose(output_file);

    return 0;
}