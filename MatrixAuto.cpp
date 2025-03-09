#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <cstdlib>
#include <papi.h>
#include <vector>
#include <fstream>
#include <unistd.h>

using namespace std;

#define SYSTEMTIME clock_t

void saveResult(int type, int size, int cores, string time, long long values[])
{
    std::string line = to_string(type) + ',' +
                       to_string(size) + ',' +
                       to_string(cores) + ',' +
                       time;

    for (int i = 0; i < 5; i++)
    {
        line += ',' + std::to_string(values[i]); // Dereference pointer to get value
    }

    std::ofstream file("results.txt", std::ios::app);

    if (file.is_open())
    {
        file << line << "\n";
        file.close();
        cout << line << endl;
    }
    else
    {
        cout << "Error: " << line << endl;
    }
}

string OnMult(int m_ar, int m_br, int cores)
{
    double Time1, Time2;

    char st[100];
    double temp;
    int i, j, k;

    double *pha, *phb, *phc;

    pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

    for (i = 0; i < m_ar; i++)
        for (j = 0; j < m_ar; j++)
            pha[i * m_ar + j] = (double)1.0;

    for (i = 0; i < m_br; i++)
        for (j = 0; j < m_br; j++)
            phb[i * m_br + j] = (double)(i + 1);

    omp_set_num_threads(cores);

    Time1 = omp_get_wtime();

#pragma omp parallel for private(j, k, temp) schedule(static)
    for (i = 0; i < m_ar; i++)
    {
        for (j = 0; j < m_br; j++)
        {
            temp = 0.0;
            for (k = 0; k < m_ar; k++)
            {
                temp += pha[i * m_ar + k] * phb[k * m_br + j];
            }
            phc[i * m_ar + j] = temp;
        }
    }

    Time2 = omp_get_wtime();
    sprintf(st, "Time: %.3f seconds\n", (double)(Time2 - Time1));
    cout << st;
    sprintf(st, "%.3f", (double)(Time2 - Time1));

    cout << "Result matrix: " << endl;
    for (i = 0; i < 1; i++)
    {
        for (j = 0; j < min(10, m_br); j++)
            cout << phc[j] << " ";
    }
    cout << endl;

    free(pha);
    free(phb);
    free(phc);

    return st;
}

string OnMultLine(int m_ar, int m_br, int cores)
{
    double Time1, Time2;

    char st[100];
    double temp;
    int i, j, k;

    double *pha, *phb, *phc;

    pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
    phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

    for (i = 0; i < m_ar; i++)
        for (j = 0; j < m_ar; j++)
            pha[i * m_ar + j] = (double)1.0;

    for (i = 0; i < m_br; i++)
        for (j = 0; j < m_br; j++)
            phb[i * m_br + j] = (double)(i + 1);

    omp_set_num_threads(cores);

    Time1 = omp_get_wtime();

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < m_ar; i++)
    {
        for (int j = 0; j < m_br; j++) 
        {
            for (int k = 0; k < m_br; k++)
            {
                phc[i * m_br + k] += pha[i * m_br + j] * phb[j * m_br + k];
            }
        }
    }

    Time2 = omp_get_wtime();

    std::cout << std::fixed << std::setprecision(3) << "Time: " << (Time2 - Time1) << " seconds\n";
    sprintf(st, "Time: %3.3f seconds\n", (double)(Time2 - Time1));
    cout << st;
    sprintf(st, "%3.3f", (double)(Time2 - Time1));

    cout << "Result matrix: " << endl;
    for (i = 0; i < 1; i++)
    {
        for (j = 0; j < min(10, m_br); j++)
            cout << phc[j] << " ";
    }
    cout << endl;

    free(pha);
    free(phb);
    free(phc);

    return st;
}

float produtoInterno(float *v1, float *v2, int col)
{
    int i;
    float soma = 0.0;

    for (i = 0; i < col; i++)
        soma += v1[i] * v2[i];

    return (soma);
}

void handle_error(int retval)
{
    printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
    exit(1);
}

void init_papi()
{
    int retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT && retval < 0)
    {
        printf("PAPI library version mismatch!\n");
        exit(1);
    }
    if (retval < 0)
        handle_error(retval);

    std::cout << "PAPI Version Number: MAJOR: " << PAPI_VERSION_MAJOR(retval)
              << " MINOR: " << PAPI_VERSION_MINOR(retval)
              << " REVISION: " << PAPI_VERSION_REVISION(retval) << "\n";
}

int caller(int type, int size, int cores)
{
    int EventSet = PAPI_NULL;
    long long values[5];
    int ret;

    ret = PAPI_library_init(PAPI_VER_CURRENT);
    if (ret != PAPI_VER_CURRENT)
        std::cout << "FAIL" << endl;

    ret = PAPI_create_eventset(&EventSet);
    if (ret != PAPI_OK)
        cout << "ERRO: create eventset" << endl;

    ret = PAPI_add_event(EventSet, PAPI_L1_DCM);
    if (ret != PAPI_OK)
        cout << "ERRO: PAPI_L1_DCM" << endl;

    ret = PAPI_add_event(EventSet, PAPI_L2_DCM);
    if (ret != PAPI_OK)
        cout << "ERRO: PAPI_L2_DCM" << endl;

    ret = PAPI_add_event(EventSet, PAPI_FP_INS); // Floating point instructions
    if (ret != PAPI_OK)
        cout << "ERROR: PAPI_FP_INS" << endl;

    ret = PAPI_add_event(EventSet, PAPI_TOT_INS); // Total instructions completed
    if (ret != PAPI_OK)
        cout << "ERROR: PAPI_TOT_INS" << endl;

    ret = PAPI_add_event(EventSet, PAPI_FP_OPS); // Total cycles completed
    if (ret != PAPI_OK)
        cout << "ERROR: PAPI_FP_OPS" << endl;

    ret = PAPI_start(EventSet);
    if (ret != PAPI_OK)
        cout << "ERRO: PAPI_START" << endl;

    string time;

    switch (type)
    {
    case 1:
        time = OnMult(size, size, cores);
        break;
    case 2:
        time = OnMultLine(size, size, cores);
        break;
    default:
        break;
    }

    ret = PAPI_stop(EventSet, values);
    if (ret != PAPI_OK)
        cout << "ERRO: Stop PAPI" << endl;

    saveResult(type, size, cores, time, values);

    ret = PAPI_reset(EventSet);
    if (ret != PAPI_OK)
        std::cout << "FAIL reset" << endl;

    ret = PAPI_remove_event(EventSet, PAPI_L1_DCM);
    if (ret != PAPI_OK)
        std::cout << "FAIL remove event" << endl;

    ret = PAPI_remove_event(EventSet, PAPI_L2_DCM);
    if (ret != PAPI_OK)
        std::cout << "FAIL remove event" << endl;

    ret = PAPI_remove_event(EventSet, PAPI_FP_INS); // Floating point instructions
    if (ret != PAPI_OK)
        cout << "ERROR: remove PAPI_FP_INS" << endl;

    ret = PAPI_remove_event(EventSet, PAPI_TOT_INS); // Total instructions completed
    if (ret != PAPI_OK)
        cout << "ERROR: remove PAPI_TOT_INS" << endl;

    ret = PAPI_remove_event(EventSet, PAPI_FP_OPS); // Total cycles completed
    if (ret != PAPI_OK)
        cout << "ERROR: remove PAPI_FP_OPS" << endl;

    ret = PAPI_destroy_eventset(&EventSet);
    if (ret != PAPI_OK)
        std::cout << "FAIL destroy" << endl;

    return 0;
}

int main(int argc, char *argv[])
{
    vector<int> size = {};
    int i = 200;
    while (i < 3000)
    {
        i += 400;
        size.push_back(i);
    }

    vector<int> size2 = {};
    i = 4096 - 2048;
    while (i < 10240)
    {
        i += 2048;
        size2.push_back(i);
    }

    vector<int> cores = {};
    i = 1;
    while (i < 16)
    {
        i++;
        cores.push_back(i);
    }

    for (int i : size2) {
        caller(2, i, 10);
    }

    for (int i : size2) {
        caller(2, i, 16);
    }
    
    /*
    for (int i : size2)
    {
        caller(2, i, 1);
    }

    
    for (int c : cores) {
        for (int i : size) {
            caller(1, i, c);
        }
    }

    for (int c : cores) {
        for (int i : size) {
            caller(2, i, c);
        }
    }


    for (int i : size)
    {
        caller(1, i, 1);
    }

    for (int i : size)
    {
        caller(2, i, 1);
    }

    for (int i : size2)
    {
        caller(2, i, 1);
    }

    while (true)
    {
        cout << '\a';
    }*/

    return 0;
}
