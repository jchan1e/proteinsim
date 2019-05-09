
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include <iostream>
#include <fstream>
#include <vector>

//#include "CSCIx229.h"
//#include <SDL.h>
//#include <SDL_opengl.h>
//#include "objects.h"

using namespace std;

////////////////////

int compMode = 0; // 0 = single thread
                  // 1 = multi thread
                  // 2 = GPU

///////////////////////////////////

__global__ void physics(const int n, const int frame, const float* aminosPrev, float* aminosNext, float* history)
{
  float k = 1.0; // Bond Spring Constant
  float ke = -0.01; // Electrostatic Constant
  float kh = -0.2; // Hydrophobicity Constant
  float kc = 1.0;  // Collision Force Constant

  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n)
  {
    float hx = 0.0;
    float hy = 0.0;
    float hz = 0.0;
    float ex = 0.0;
    float ey = 0.0;
    float ez = 0.0;
    float cx = 0.0;
    float cy = 0.0;
    float cz = 0.0;
    float fx = 0.0;
    float fy = 0.0;
    float fz = 0.0;
    for (int j=0; j < n; ++j)
    {
      if (i != j)
      {
        //calculate forces
        float x1 = aminosPrev[8*i + 0];
        float y1 = aminosPrev[8*i + 1];
        float z1 = aminosPrev[8*i + 2];
        float x2 = aminosPrev[8*j + 0];
        float y2 = aminosPrev[8*j + 1];
        float z2 = aminosPrev[8*j + 2];
        float dx = x2-x1;
        float dy = y2-y1;
        float dz = z2-z1;
        float dist = sqrt(dx*dx + dy*dy + dz*dz);
        float h1 = aminosPrev[8*i + 6];
        float e1 = aminosPrev[8*i + 7];
        float h2 = aminosPrev[8*j + 6];
        float e2 = aminosPrev[8*j + 7];

        float vx = dx/dist;
        float vy = dy/dist;
        float vz = dz/dist;

        // Hydrophobic forces
        // Fh = Kh*h1*h2/(r^14-r^8)
        float d = max(dist, 1.0f);
        hx = kh*h1*h2*(pow(d,-14) - pow(d,-8)) * vx;
        hy = kh*h1*h2*(pow(d,-14) - pow(d,-8)) * vy;
        hz = kh*h1*h2*(pow(d,-14) - pow(d,-8)) * vz;

        // Electrsostatic forces
        // Fe = k*q1*q2/r^2
        ex = ke*e1*e2/min(dist*dist, 1.0f) * vx;
        ey = ke*e1*e2/min(dist*dist, 1.0f) * vy;
        ez = ke*e1*e2/min(dist*dist, 1.0f) * vz;

        // Collision forces
        // soft collisions, spring force model
        cx = 0.0;
        cy = 0.0;
        cz = 0.0;
        if (dist < 1.0)
        {
          cx = kc*(1.0-dist) * vx;
          cy = kc*(1.0-dist) * vy;
          cz = kc*(1.0-dist) * vz;
        }

        //if (hx*hx + hy*hy + hz*hz > 0.01)
        //  cout << "H" << i << ":\t"<< hx << "\t" << hy << "\t" << hz << "\t" << dist << endl;
        //if (ex*ex + ey*ey + ez*ez > 0.01)
        //  cout << "E" << i << ":\t"<< ex << "\t" << ey << "\t" << ez << "\t" << dist << endl;
        //if (cx*cx + cy*cy + cz*cz > 0.01)
        //  cout << "C" << i << ":\t"<< cx << "\t" << cy << "\t" << cz << "\t" << dist << endl;

        fx += hx + ex - cx;
        fy += hy + ey - cy;
        fz += hz + ez - cz;
      }
    }
    // spring tension
    if (i > 0)
    {
      float x1 = aminosPrev[8*i + 0];
      float y1 = aminosPrev[8*i + 1];
      float z1 = aminosPrev[8*i + 2];
      float x2 = aminosPrev[8*(i-1) + 0];
      float y2 = aminosPrev[8*(i-1) + 1];
      float z2 = aminosPrev[8*(i-1) + 2];
      float dx = x2-x1;
      float dy = y2-y1;
      float dz = z2-z1;
      float dist = sqrt(dx*dx + dy*dy + dz*dz);
      fx += k*(dist-1.0) * dx/dist;
      fy += k*(dist-1.0) * dy/dist;
      fz += k*(dist-1.0) * dz/dist;
    }
    if (i < n-1)
    {
      float x1 = aminosPrev[8*i + 0];
      float y1 = aminosPrev[8*i + 1];
      float z1 = aminosPrev[8*i + 2];
      float x2 = aminosPrev[8*(i+1) + 0];
      float y2 = aminosPrev[8*(i+1) + 1];
      float z2 = aminosPrev[8*(i+1) + 2];
      float dx = x2-x1;
      float dy = y2-y1;
      float dz = z2-z1;
      float dist = sqrt(dx*dx + dy*dy + dz*dz);
      fx += k*(dist-1.0) * dx/dist;
      fy += k*(dist-1.0) * dy/dist;
      fz += k*(dist-1.0) * dz/dist;
    }

    // update velocities
    aminosNext[8*i+3] = aminosPrev[8*i+3] + fx;
    aminosNext[8*i+4] = aminosPrev[8*i+4] + fy;
    aminosNext[8*i+5] = aminosPrev[8*i+5] + fz;
    // damping
    aminosNext[8*i+3] *= 0.9995;
    aminosNext[8*i+4] *= 0.9995;
    aminosNext[8*i+5] *= 0.9995;
    // update positions
    aminosNext[8*i+0] = aminosPrev[8*i+0] + 0.1*aminosNext[8*i+3];
    aminosNext[8*i+1] = aminosPrev[8*i+1] + 0.1*aminosNext[8*i+4];
    aminosNext[8*i+2] = aminosPrev[8*i+2] + 0.1*aminosNext[8*i+5];
    // copy to history
    history[3*n*frame + 3*i + 0] = aminosNext[8*i+0];
    history[3*n*frame + 3*i + 1] = aminosNext[8*i+1];
    history[3*n*frame + 3*i + 2] = aminosNext[8*i+2];
  }
}

///////////////////////////////////

//void dummy(int n, float* nodes, float* hist)
//{
//  for (int i=0; i < n; ++i)
//  {
//    float hx = 0.0;
//    float hy = 0.0;
//    float hz = 0.0;
//    float ex = 0.0;
//    float ey = 0.0;
//    float ez = 0.0;
//    float cx = 0.0;
//    float cy = 0.0;
//    float cz = 0.0;
//    float fx = 0.0;
//    float fy = 0.0;
//    float fz = 0.0;
//    for (int j=0; j < n; ++j)
//    {
//      if (i != j)
//      {
//        float x1 = nodes[8*i + 0];
//        float y1 = nodes[8*i + 1];
//        float z1 = nodes[8*i + 2];
//        float x2 = nodes[8*j + 0];
//        float y2 = nodes[8*j + 1];
//        float z2 = nodes[8*j + 2];
//        float dx = x2-x1;
//        float dy = y2-y1;
//        float dz = z2-z1;
//        float dist = sqrt(dx*dx + dy*dy + dz*dz);
//        float h1 = nodes[8*i + 6];
//        float e1 = nodes[8*i + 7];
//        float h2 = nodes[8*j + 6];
//        float e2 = nodes[8*j + 7];
//
//        // Hydrophobic forces
//        // Fh = Kh*h1*h2/(r^14-r^8)
//        float d = max(dist, 1.0f);
//        hx = kh*h1*h2*(pow(d,-14) - pow(d,-8)) * dx/dist;
//        hy = kh*h1*h2*(pow(d,-14) - pow(d,-8)) * dy/dist;
//        hz = kh*h1*h2*(pow(d,-14) - pow(d,-8)) * dz/dist;
//
//        // Electrsostatic forces
//        // Fe = k*q1*q2/r^2
//        ex = ke*e1*e2/min(dist*dist, 1.0f) * dx/dist;
//        ey = ke*e1*e2/min(dist*dist, 1.0f) * dy/dist;
//        ez = ke*e1*e2/min(dist*dist, 1.0f) * dz/dist;
//
//        // Collision forces
//        // soft collisions, spring force model
//        if (dist < 1.0)
//        {
//          cx = kc*(1.0-dist) * dx/dist;
//          cy = kc*(1.0-dist) * dy/dist;
//          cz = kc*(1.0-dist) * dz/dist;
//        }
//        else
//        {
//          cx = 0.0;
//          cy = 0.0;
//          cz = 0.0;
//        }
//
//        //if (hx*hx + hy*hy + hz*hz > 0.01)
//        //  cout << "H" << i << ":\t"<< hx << "\t" << hy << "\t" << hz << "\t" << dist << endl;
//        //if (ex*ex + ey*ey + ez*ez > 0.01)
//        //  cout << "E" << i << ":\t"<< ex << "\t" << ey << "\t" << ez << "\t" << dist << endl;
//        //if (cx*cx + cy*cy + cz*cz > 0.01)
//        //  cout << "C" << i << ":\t"<< cx << "\t" << cy << "\t" << cz << "\t" << dist << endl;
//
//        fx += hx + ex - cx;
//        fy += hy + ey - cy;
//        fz += hz + ez - cz;
//      }
//    }
//    // update velocities
//    nodes[8*i + 3] += fx;
//    nodes[8*i + 4] += fy;
//    nodes[8*i + 5] += fz;
//  }
//  // Spring Tension
//  for (int i=0; i < n-1; ++i)
//  {
//    int j = i + 1;
//    float x1 = nodes[8*i + 0];
//    float y1 = nodes[8*i + 1];
//    float z1 = nodes[8*i + 2];
//    float x2 = nodes[8*j + 0];
//    float y2 = nodes[8*j + 1];
//    float z2 = nodes[8*j + 2];
//    float dx = x2-x1;
//    float dy = y2-y1;
//    float dz = z2-z1;
//    float dist = sqrt(dx*dx + dy*dy + dz*dz);
//    nodes[8*i + 3] += k*(dist-1.0) * dx/dist;
//    nodes[8*i + 4] += k*(dist-1.0) * dy/dist;
//    nodes[8*i + 5] += k*(dist-1.0) * dz/dist;
//    nodes[8*j + 3] -= k*(dist-1.0) * dx/dist;
//    nodes[8*j + 4] -= k*(dist-1.0) * dy/dist;
//    nodes[8*j + 5] -= k*(dist-1.0) * dz/dist;
//
//    //if (dist < 0.9 || dist > 1.25)
//    //  cout << dist << endl;
//  }
//  for (int i=0; i < n; ++i)
//  {
//    // damping
//    nodes[8*i + 3] *= 0.9995;
//    nodes[8*i + 4] *= 0.9995;
//    nodes[8*i + 5] *= 0.9995;
//    // update positions
//    nodes[8*i + 0] += 0.1*nodes[8*i + 3];
//    nodes[8*i + 1] += 0.1*nodes[8*i + 4];
//    nodes[8*i + 2] += 0.1*nodes[8*i + 5];
//    //hist->push_back(nodes[8*i + 0]);
//    //hist->push_back(nodes[8*i + 1]);
//    //hist->push_back(nodes[8*i + 2]);
//  }
//}

int main(int argc, char *argv[])
{
  // flags
  for (int i=0; i < argc; ++i)
  {
    if (strcmp(argv[i],"-m") == 0)
    {
      compMode = 1;
      for (int j=i+1; j < argc; ++j)
      {
        argv[j-1] = argv[j];
      }
      i--;
      argc--;
    }
    if (strcmp(argv[i],"-g") == 0)
    {
      compMode = 2;
      for (int j=i+1; j < argc; ++j)
      {
        argv[j-1] = argv[j];
      }
      i--;
      argc--;
    }
  }

  // args
  if (argc != 3 && argc != 4)
  {
    cerr << "Usage: sim infile outfile [num_frames]\n";
    return 1;
  }

  //Initialize
  int num_frames = 6000;
  if (argc == 4) num_frames = stoi(argv[3]);
  ifstream infile(argv[1]);
  if (!infile.is_open())
  {
    cerr << "could not open file " << argv[1] << endl;
    return 1;
  }
  string line;
  getline(infile, line);
  int nAminos = stoi(line);
  cout << nAminos << endl;
  float* aminos = new float[nAminos*8];
  for (int i=0; i < nAminos; ++i)
  {
    aminos[8*i + 0] = nAminos/2.0 - i; // x coordinate
    aminos[8*i + 1] = 0.0;           // y coordinate
    aminos[8*i + 2] = 0.0;           // z coordinate
    aminos[8*i + 3] = 0.0;           // x velocity
    aminos[8*i + 4] = 0.0;           // y velocity
    aminos[8*i + 5] = 0.0;           // z velocity
    getline(infile, line);
    aminos[8*i + 6] = stof(line);    // hydrophobicity
    getline(infile, line);
    aminos[8*i + 7] = stof(line);    // electrostatic charge
  }
  infile.close();
  aminos[1] = 0.1;
  aminos[5] = 0.01;
  //aminos[nAminos*8-7] = -0.1;
  //aminos[nAminos*8-3] = -0.01;

  //vector<float> history;
  float* history = new float[3*nAminos*num_frames];
  int frames = 0;

  float* g_aminos1 = NULL;
  float* g_aminos2 = NULL;
  float* g_history = NULL;
  cudaMalloc(&g_aminos1, 8*nAminos*sizeof(float));
  cudaMalloc(&g_aminos2, 8*nAminos*sizeof(float));
  cudaMalloc(&g_history, 3*num_frames*nAminos*sizeof(float));

  cudaMemcpy(g_aminos1, aminos, 8*nAminos*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(g_aminos2, aminos, 8*nAminos*sizeof(float), cudaMemcpyHostToDevice);

  ////////Main Loop////////
  while (frames < num_frames)
  {
    physics<<<max(1,nAminos/128),128>>>(nAminos, frames, g_aminos1, g_aminos2, g_history);
//    float test[8];
//    cudaMemcpy(&test, g_aminos2, 8*sizeof(float), cudaMemcpyDeviceToHost);
//    for (int i=0; i < 8; ++i)
//      cout << test[i] << " ";
//    cout << endl;
    frames += 1;
    float* tmp = g_aminos1;
    g_aminos1 = g_aminos2;
    g_aminos2 = tmp;
  }

  // wait and copy back history
  cudaMemcpy(history, g_history, 3*num_frames*nAminos*sizeof(float), cudaMemcpyDeviceToHost);

//  for (int i=0; i < frames; ++i)
//  {
//    cout << history[3*nAminos*i + 0 + 0] << " ";
//    cout << history[3*nAminos*i + 0 + 1] << " ";
//    cout << history[3*nAminos*i + 0 + 2] << endl;
//  }

  // write to file
  ofstream outfile(argv[2], ofstream::binary);
  if (!outfile.is_open())
    cerr << "could not open file: " << argv[2] << endl;
  else
  {
    float aminoList[2*nAminos];
    for (int i=0; i < nAminos; ++i)
    {
      aminoList[2*i]   = aminos[8*i + 6];
      aminoList[2*i+1] = aminos[8*i + 7];
    }
    outfile.write((char*)&nAminos, sizeof(int));
    outfile.write((char*)&frames, sizeof(int));
    outfile.write((char*)aminoList, 2*nAminos*sizeof(float));
    outfile.write((char*)history, 3*frames*nAminos*sizeof(float));
    outfile.close();
  }

  //cout << "Shutting Down\n";

  delete[] aminos;
  delete[] history;
  cudaFree(g_aminos1);
  cudaFree(g_aminos2);
  cudaFree(g_history);

  return 0;
}
