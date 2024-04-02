#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include "melpe.h"
#include "sc1200.h"
#include "npp.h"

#define FRAME_SIZE 180 // 将每帧的样本数更改为 540
#define BITS_PER_FRAME 54 // 每帧的比特数，根据帧长计算
#define SAMPLE_RATE 8000 // 采样率为 8kHz
#define BITS_PER_SAMPLE 16

struct melp_param pa;
struct melp_param pat;

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

int processWavFile(const char *inputPath, const char *outputPath) {
    // 在这里实现 WAV 文件处理的逻辑，根据 WAV 文件格式修

    FILE *input_file, *output_file;
    short speech[FRAME_SIZE];
    unsigned char bitstream[(BITS_PER_FRAME + 7) / 8];
    WAVHeader header;

    input_file = fopen(inputPath, "rb");
    if (input_file == NULL) {
        printf("无法打开输入文件：%s\n", inputPath);
        return 1;
    }

    output_file = fopen(outputPath, "wb");
    if (output_file == NULL) {
        printf("无法打开输出文件：%s\n", outputPath);
        fclose(input_file);
        return 1;
    }

    melpe_i();

    fseek(output_file, sizeof(WAVHeader), SEEK_SET);
    uint32_t dataSize = 0;

    while (fread(speech, sizeof(short), FRAME_SIZE, input_file) == FRAME_SIZE) {
        analysis(speech, &pa);
        rmquan(&pat, speech, &pa);
        fwrite(speech, sizeof(short), FRAME_SIZE, output_file);
        dataSize += FRAME_SIZE * sizeof(short);
    }

    fseek(output_file, 0, SEEK_SET);
    createWAVHeader(&header, dataSize);
    fwrite(&header, sizeof(WAVHeader), 1, output_file);


    fclose(input_file);
    fclose(output_file);
}

void processWavFilesInDirectory(const char *inputDir, const char *outputDir) {
    DIR *dir;
    struct dirent *entry;

    dir = opendir(inputDir);

    if (dir == NULL) {
        perror("Error opening directory"); // 检查目录打开是否正常
        exit(EXIT_FAILURE);
    }

    while ((entry = readdir(dir)) != NULL) {
        if (entry->d_type == DT_REG) { // 检查读取到的是否为普通文件，DT_REG表示该文件是普通文件
            char inputPath[512]; // 存储输入文件的路径
            char outputPath[512]; // 存储输出文件的路径

            snprintf(inputPath, sizeof(inputPath), "%s/%s", inputDir, entry->d_name);
            inputPath[sizeof(inputPath) - 1] = '\0';  // 确保字符串以 null 结尾
            snprintf(outputPath, sizeof(outputPath), "%s/%s", outputDir, entry->d_name);
            outputPath[sizeof(outputPath) - 1] = '\0';  // 确保字符串以 null 结尾


            processWavFile(inputPath, outputPath);
       }
    }

    closedir(dir);
}

int main() {
    const char *inputDirectory = "/mnt/e/data/ls8k"; // 输入文件夹路径
    const char *outputDirectory = "/mnt/e/data/ls_2400_rq_8k"; // 输出文件夹路径

    processWavFilesInDirectory(inputDirectory, outputDirectory);
    return 0;
}
