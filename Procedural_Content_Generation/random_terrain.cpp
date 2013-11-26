#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>

#define max_dim 32769*32769

using namespace std;

class NODE{
  public:
    NODE(){return;}
    NODE(unsigned int, unsigned int);
    ~NODE(){}
    unsigned int xpos;
    unsigned int ypos;
    NODE *parent;
} *trunk, *cur_node, *dummy;

NODE::NODE(unsigned int new_x,unsigned int new_y)
{
  xpos = new_x;
  ypos = new_y;
}

unsigned char*** RGB_heights;
unsigned char** height_map;
unsigned int width,height,image_size,file_size,buffer_bt;
bool** Done_Yet;
ofstream outfile;

void write_header();

float frand(){return (rand()/(float)RAND_MAX);}

void full_noise();
void diamond_square_color();
void float_ds_color(float);
bool diamond_square();
bool diamond_square_wrap();
bool float_ds(float);
void water_fill(int, int, unsigned char);
void water_fill(int, int, unsigned char,unsigned long int);
bool water_bodies(int);
void basic_avg(float);
void adv_avg(float);
void hght_fill(float);
void seeded_dsa();
void TwoD_MD();

int main(int argc, char* argv[])
{
  srand(time(NULL));
  char filename[24];
  unsigned char method;
  float scale = .7;
  if (argc<2)
  {
    cout<<"Enter width then height separated by space. Ex: '255 300'. ";
    cin>>width>>height;
    cout<<"Enter filename including '.bmp'. Ex: 'picture.bmp'. ";
    cin>>filename;
    cout<<"Enter single-character identifier for method desired then 0 if not doing a float method."<<endl;
    cout<<"Currently implemented:"<<endl;
    cout<<"\tb\tdiamond-square method with wrapping"<<endl;
    cout<<"\tc\tdiamond-square colors"<<endl;
    cout<<"\td\tdiamond-square method"<<endl;
    cout<<"\te\tdiamond-square colors by floats (include constant between 0 and 1)"<<endl;
    cout<<"\tf\tdiamond-square method by floats (include constant between 0 and 1)"<<endl;
    cout<<"\tg\tseeded dsa with color"<<endl;
    cout<<"\ts\trandom noise method"<<endl;
    cout<<"\tv\tmidpoint displacement"<<endl;
    cout<<"\tw\twater by seeds (include integer)"<<endl;
    cout<<"\tx\twater by height"<<endl;
    cout<<"\ty\tbasic averaging"<<endl;
    cout<<"\tz\tadvanced averaging"<<endl;
    cin>>method>>scale;
  }
  else if (argc==5)
  {
    width = atoi(argv[1]);
    height = atoi (argv[2]);
    sprintf(filename,"%s",argv[3]);
    method = argv[4][0];
  }
  unsigned int num_cells = width*height;
  buffer_bt = (((width*24+31)/32) * 4) - (width*3);
  if (num_cells > max_dim) {cout<<"Image size too big! Must be less than "<<max_dim<<"."<<endl; return 0;}
  image_size = num_cells * 3 + buffer_bt * height;
  file_size = image_size + 54;
  height_map = new unsigned char*[height];
  for (unsigned int i = 0; i < height; i++) height_map[i] = new unsigned char[width];
  for (unsigned int i = 0; i < height; i++)
  {
    for (unsigned int j = 0; j < width; j++) height_map[i][j] = 128;
  }
  bool ran_well = true;
  bool done = false;
  switch (method)
  {
    case 'c':
      outfile.open((const char*)filename);
      write_header();
      diamond_square_color();
      outfile.close();
      done = true;
      break;
    case 'e':
      outfile.open((const char*)filename);
      write_header();
      float_ds_color(scale);
      outfile.close();
      done = true;
      break;
    case 'g':
      outfile.open((const char*)filename);
      write_header();
      seeded_dsa();
      outfile.close();
      done = true;
      break;
    case 'v':
      outfile.open((const char*)filename);
      write_header();
      TwoD_MD();
      outfile.close();
      done = true;
      break;
    case 'w':
      outfile.open((const char*)filename);
      write_header();
      water_bodies((int)scale+2);
      outfile.close();
      done = true;
      break;
    case 'x':
      outfile.open((const char*)filename);
      write_header();
      hght_fill(scale);
      outfile.close();
      done = true;
      break;
    case 'y':
      outfile.open((const char*)filename);
      write_header();
      basic_avg(scale);
      outfile.close();
      done = true;
      break;
    case 'z':
      outfile.open((const char*)filename);
      write_header();
      adv_avg(scale);
      outfile.close();
      done = true;
      break;
    default:
      break;
  }
  if (done) return 0;
  switch(method)
  {
    case 'b':
      ran_well = diamond_square_wrap();
      break;
    case 'd':
      ran_well = diamond_square();
      break;
    case 'f':
      ran_well = float_ds(scale);
      break;
    case 's':
      full_noise();
      break;
    default:
      ran_well = false;
      break;
  }
  if (ran_well)
  {
    outfile.open((const char*)filename);
    write_header();
    for (unsigned int i = 0; i < height; i++)
    {
      for (unsigned int j = 0; j < width; j++) outfile << height_map[i][j] << height_map[i][j] << height_map[i][j];
      for (unsigned int j = 0; j < buffer_bt; j++) outfile << (char)0;
      delete height_map[i];
    }
    delete height_map;
    outfile.close();
  }
  else cout<<"Failure! See messages."<<endl;
  return 0;
}

bool diamond_square()
{
  if ((width & 1)==0) {cout<<"Bad width for DSA! Must be odd."; return 0;}
  else if ((height & 1)==0) {cout<<"Bad height for DSA! Must be odd."; return 0;}
  else
  {
    for (unsigned int i = 1; i < 32; i++)
    {
      if ((1<<(i-1))==(width>>1)) break;
      else if (i==31) {cout<<"Bad width for DSA! Must be 1 plus a power of 2."; return 0;}
    }
    for (unsigned int i = 1; i < 32; i++)
    {
      if ((1<<(i-1))==(height>>1)) break;
      else if (i==31) {cout<<"Bad height for DSA! Must be 1 plus a power of 2."; return 0;}
    }
    if (height != width) {cout<<"Bad dimensions for DSA! Must have w=h."; return 0;}
  }
  unsigned int sq_width = width - 1;
  unsigned int max_rand = 128;
  height_map[0][0] = height_map[width-1][0] = height_map[0][width-1] = height_map[width-1][width-1] = 128;
  int x_pos = 0;
  int y_pos = 0;
  int temp;
  int temp2;
  while (sq_width>0)
  {
    for (unsigned int i = 0; i < (width - 1); i += sq_width)
    {
      for (unsigned int j = 0; j < (width - 1); j += sq_width)
      {
        x_pos = (i + (sq_width/2));
        y_pos = (j + (sq_width/2));
        height_map[y_pos][x_pos] = (unsigned char)max((unsigned)0,min((unsigned)255,(height_map[j][i] + height_map[j+sq_width][i] + height_map[j][i+sq_width] + height_map[j+sq_width][i+sq_width])/4 + ((rand()%max_rand - (max_rand/2)))));
      }
    }
    int k,l,m;
    m = sq_width/2;
    for (unsigned int i = 0; i < (width - 1); i += sq_width)
    {
      for (unsigned int j = 0; j < (width - 1); j += sq_width)
      {
        x_pos = (i + m);
        y_pos = (j + m);
        //do up first
        l = x_pos;
        k = j;
//        height_map[k][l] = (unsigned char)max((unsigned int)0,min((unsigned int)255,(unsigned int)((int)((height_map[(k-m<0?k-m+width:k-m)][l] + height_map[(k+m>width?k+m-width:k+m)][l] + height_map[k][(l+m>width?l+m-width:l+m)] + height_map[k][(l-m<0?l-m+width:l-m)])/4) + ((rand()%max_rand - (max_rand/2))))));
        temp2 = rand()%max_rand - max_rand/2;
        temp = (signed int)(height_map[(k-m<0?k-m+width:k-m)][l] + height_map[(k+m>width?k+m-width:k+m)][l] + height_map[k][(l+m>width?l+m-width:l+m)] + height_map[k][(l-m<0?l-m+width:l-m)])/4 + temp2;
        if (temp>255) height_map[k][l] = 255;
        else if (temp<0) height_map[k][l] = 0;
        else height_map[k][l] = temp;
        //do left next
        l = i;
        k = y_pos;
        temp2 = rand()%max_rand - max_rand/2;
        temp = (signed int)(height_map[(k-m<0?k-m+width:k-m)][l] + height_map[(k+m>width?k+m-width:k+m)][l] + height_map[k][(l+m>width?l+m-width:l+m)] + height_map[k][(l-m<0?l-m+width:l-m)])/4 + temp2;
        if (temp>255) height_map[k][l] = 255;
        else if (temp<0) height_map[k][l] = 0;
        else height_map[k][l] = temp;
        //check and do right then down
        if ((x_pos + sq_width)>width)
        {
          l = i+sq_width;
          k = y_pos;
          temp2 = rand()%max_rand - max_rand/2;
          temp = (signed int)(height_map[(k-m<0?k-m+width:k-m)][l] + height_map[(k+m>width?k+m-width:k+m)][l] + height_map[k][(l+m>width?l+m-width:l+m)] + height_map[k][(l-m<0?l-m+width:l-m)])/4 + temp2;
          if (temp>255) height_map[k][l] = 255;
          else if (temp<0) height_map[k][l] = 0;
          else height_map[k][l] = temp;
        }
        if ((y_pos + sq_width)>width)
        {
          l = x_pos;
          k = j + sq_width;
          temp2 = rand()%max_rand - max_rand/2;
          temp = (signed int)(height_map[(k-m<0?k-m+width:k-m)][l] + height_map[(k+m>width?k+m-width:k+m)][l] + height_map[k][(l+m>width?l+m-width:l+m)] + height_map[k][(l-m<0?l-m+width:l-m)])/4 + temp2;
          if (temp>255) height_map[k][l] = 255;
          else if (temp<0) height_map[k][l] = 0;
          else height_map[k][l] = temp;
        }
      }
    }
    sq_width /= 2;
    max_rand = (max_rand>2?max_rand/2:2);
  }
  return 1;
}

void full_noise()
{
  int k=8;
  while (k--)
  {
    for (unsigned int i = 0; i < height; i++)
    {
      for (unsigned int j = 0; j < width; j++) height_map[i][j] += (char)((rand() % (2<<k)) - (2<<k));
    }
  }
}

void write_header()
{
  unsigned char place_hold;
  outfile << "B" << "M";
  for (unsigned int i = 0; i < 32; i += 8)
  {
    outfile << (char)((file_size >> i)%256);
  }
  place_hold = 0;
  outfile << place_hold << place_hold << place_hold << place_hold;
  place_hold = 54;
  outfile << place_hold;
  place_hold = 0;
  outfile << place_hold << place_hold << place_hold;
  place_hold = 40;
  outfile << place_hold;
  place_hold = 0;
  outfile << place_hold << place_hold << place_hold;
  for (unsigned int i = 0; i < 32; i += 8)
  {
    outfile << (char)((width >> i)%256);
  }
  for (unsigned int i = 0; i < 32; i += 8)
  {
    outfile << (char)((height >> i)%256);
  }
  place_hold = 1;
  outfile << place_hold;
  place_hold = 0;
  outfile << place_hold;
  place_hold = 24;
  outfile << place_hold;
  place_hold = 0;
  for (unsigned int i = 0; i < 5; i++) outfile << place_hold;
  for (unsigned int i = 0; i < 32; i+= 8)
  {
    outfile << (char)((image_size >> i)%256);
  }
  for (unsigned int i = 0; i < 32; i+= 8)
  {
    outfile << (char)((3780 >> i)%256);
  }
  for (unsigned int i = 0; i < 32; i+= 8)
  {
    outfile << (char)((3780 >> i)%256);
  }
  for (unsigned int i = 0; i < 8; i++)
  {
    outfile << (char)0;
  }
}

void diamond_square_color()
{
  RGB_heights = new unsigned char**[height];
  for (unsigned int i = 0; i < height; i++)
  {
    RGB_heights[i] = new unsigned char*[width];
    for (unsigned int j = 0; j < width; j++) RGB_heights[i][j] = new unsigned char[3];
  }
  for (int i=0;i<3;i++)
  {
    diamond_square();
    for (unsigned int j=0;j<height;j++)
    {
      for (unsigned int k=0;k<width;k++)
      {
        RGB_heights[j][k][i] = height_map[j][k];
      }
    }
  }
  for (unsigned int i = 0; i < height; i++)
  {
    for (unsigned int j = 0; j < width; j++)
    {
      for (int k=0;k<3;k++) outfile << (char)RGB_heights[i][j][k];
      delete RGB_heights[i][j];
    }
    for (unsigned int j = 0; j < buffer_bt; j++) outfile << (char)0;
    delete RGB_heights[i];
  }
}

bool diamond_square_wrap()
{
  if ((width & 1)==0) {cout<<"Bad width for DSA! Must be odd."; return 0;}
  else if ((height & 1)==0) {cout<<"Bad height for DSA! Must be odd."; return 0;}
  else
  {
    for (unsigned int i = 1; i < 32; i++)
    {
      if ((1<<(i-1))==(width>>1)) break;
      else if (i==31) {cout<<"Bad width for DSA! Must be 1 plus a power of 2."; return 0;}
    }
    for (unsigned int i = 1; i < 32; i++)
    {
      if ((1<<(i-1))==(height>>1)) break;
      else if (i==31) {cout<<"Bad height for DSA! Must be 1 plus a power of 2."; return 0;}
    }
    if (height != width) {cout<<"Bad dimensions for DSA! Must have w=h."; return 0;}
  }
  unsigned int sq_width = width - 1;
  unsigned int max_rand = 128;
  height_map[0][0] = height_map[width-1][0] = height_map[0][width-1] = height_map[width-1][width-1] = 128;
  int x_pos = 0;
  int y_pos = 0;
  int temp;
  int temp2;
  while (sq_width>0)
  {
    for (unsigned int i = 0; i < (width - 1); i += sq_width)
    {
      for (unsigned int j = 0; j < (width - 1); j += sq_width)
      {
        x_pos = (i + (sq_width/2));
        y_pos = (j + (sq_width/2));
        height_map[y_pos][x_pos] = (unsigned char)max((unsigned)0,min((unsigned)255,(height_map[j][i] + height_map[j+sq_width][i] + height_map[j][i+sq_width] + height_map[j+sq_width][i+sq_width])/4 + ((rand()%max_rand - (max_rand/2)))));
      }
    }
    int k,l,m;
    m = sq_width/2;
    for (unsigned int i = 0; i < (width - 1); i += sq_width)
    {
      for (unsigned int j = 0; j < (width - 1); j += sq_width)
      {
        x_pos = (i + m);
        y_pos = (j + m);
        //do up first
        l = x_pos;
        k = j;
//        height_map[k][l] = (unsigned char)max((unsigned int)0,min((unsigned int)255,(unsigned int)((int)((height_map[(k-m<0?k-m+width:k-m)][l] + height_map[(k+m>width?k+m-width:k+m)][l] + height_map[k][(l+m>width?l+m-width:l+m)] + height_map[k][(l-m<0?l-m+width:l-m)])/4) + ((rand()%max_rand - (max_rand/2))))));
        temp2 = rand()%max_rand - max_rand/2;
        temp = (signed int)(height_map[(k-m<0?k-m+width:k-m)][l] + height_map[(k+m>width?k+m-width:k+m)][l] + height_map[k][(l+m>width?l+m-width:l+m)] + height_map[k][(l-m<0?l-m+width:l-m)])/4 + temp2;
        if (temp>255) height_map[k][l] = 255;
        else if (temp<0) height_map[k][l] = 0;
        else height_map[k][l] = temp;
        //do left next
        l = i;
        k = y_pos;
        temp2 = rand()%max_rand - max_rand/2;
        temp = (signed int)(height_map[(k-m<0?k-m+width:k-m)][l] + height_map[(k+m>width?k+m-width:k+m)][l] + height_map[k][(l+m>width?l+m-width:l+m)] + height_map[k][(l-m<0?l-m+width:l-m)])/4 + temp2;
        if (temp>255) height_map[k][l] = 255;
        else if (temp<0) height_map[k][l] = 0;
        else height_map[k][l] = temp;
        //check and do right then down
        if ((x_pos + sq_width)>width)
        {
          l = i+sq_width;
          k = y_pos;
          height_map[k][l] = height_map[k][0];
        }
        if ((y_pos + sq_width)>width)
        {
          l = x_pos;
          k = j + sq_width;
          height_map[k][l] = height_map[0][l];
        }
      }
    }
    sq_width /= 2;
    max_rand = (max_rand>2?max_rand/2:2);
  }
  return 1;
}

bool float_ds(float factor=.5)
{
  if ((width & 1)==0) {cout<<"Bad width for DSA! Must be odd."; return 0;}
  else if ((height & 1)==0) {cout<<"Bad height for DSA! Must be odd."; return 0;}
  else
  {
    for (unsigned int i = 1; i < 32; i++)
    {
      if ((1<<(i-1))==(width>>1)) break;
      else if (i==31) {cout<<"Bad width for DSA! Must be 1 plus a power of 2."; return 0;}
    }
    for (unsigned int i = 1; i < 32; i++)
    {
      if ((1<<(i-1))==(height>>1)) break;
      else if (i==31) {cout<<"Bad height for DSA! Must be 1 plus a power of 2."; return 0;}
    }
    if (height != width) {cout<<"Bad dimensions for DSA! Must have w=h."; return 0;}
  }
  if (factor>1) factor = 1;
  else if (factor<0) factor = 0;
  unsigned int sq_width = width - 1;
  float max_rand = 32;
  float REDUCE = pow(2.0,-factor);
  float** height_flt;
  height_flt = new float*[height];
  for (unsigned int i = 0; i < height; i++) height_flt[i] = new float[width];
  height_flt[0][0] = height_flt[0][width-1] = height_flt[height-1][0] = height_flt[height-1][width-1] = 0;
  int x_pos = 0;
  int y_pos = 0;
  float temp;
  float temp2;
  while(sq_width>0)
  {
    for (unsigned int i = 0; i < (width - 1); i += sq_width)
    {
      for (unsigned int j = 0; j < (width - 1); j += sq_width)
      {
        x_pos = (i + (sq_width/2));
        y_pos = (j + (sq_width/2));
        height_flt[y_pos][x_pos] = (height_flt[j][i] + height_flt[j+sq_width][i] + height_flt[j][i+sq_width] + height_flt[j+sq_width][i+sq_width])/4 + frand()*max_rand - max_rand/2;
      }
    }
    int k,l,m;
    m = sq_width/2;
    for (unsigned int i = 0; i < (width - 1); i += sq_width)
    {
      for (unsigned int j = 0; j < (width - 1); j += sq_width)
      {
        x_pos = (i + m);
        y_pos = (j + m);
        //do up first
        l = x_pos;
        k = j;
//        height_map[k][l] = (unsigned char)max((unsigned int)0,min((unsigned int)255,(unsigned int)((int)((height_map[(k-m<0?k-m+width:k-m)][l] + height_map[(k+m>width?k+m-width:k+m)][l] + height_map[k][(l+m>width?l+m-width:l+m)] + height_map[k][(l-m<0?l-m+width:l-m)])/4) + ((rand()%max_rand - (max_rand/2))))));
        temp2 = frand()*max_rand - max_rand/2;
        temp = (height_flt[(k-m<0?k-m+width:k-m)][l] + height_flt[(k+m>width?k+m-width:k+m)][l] + height_flt[k][(l+m>width?l+m-width:l+m)] + height_flt[k][(l-m<0?l-m+width:l-m)])/4;
        height_flt[k][l] = temp + temp2;
        //do left next
        l = i;
        k = y_pos;
        temp2 = frand()*max_rand - max_rand/2;
        temp = (height_flt[(k-m<0?k-m+width:k-m)][l] + height_flt[(k+m>width?k+m-width:k+m)][l] + height_flt[k][(l+m>width?l+m-width:l+m)] + height_flt[k][(l-m<0?l-m+width:l-m)])/4;
        height_flt[k][l] = temp + temp2;
        //check and do right then down
        if ((x_pos + sq_width)>width)
        {
          l = i+sq_width;
          k = y_pos;
          temp2 = frand()*max_rand - max_rand/2;
          temp = (height_flt[(k-m<0?k-m+width:k-m)][l] + height_flt[(k+m>width?k+m-width:k+m)][l] + height_flt[k][(l+m>width?l+m-width:l+m)] + height_flt[k][(l-m<0?l-m+width:l-m)])/4;
          height_flt[k][l] = temp + temp2;
        }
        if ((y_pos + sq_width)>width)
        {
          l = x_pos;
          k = j + sq_width;
          temp2 = frand()*max_rand - max_rand/2;
          temp = (height_flt[(k-m<0?k-m+width:k-m)][l] + height_flt[(k+m>width?k+m-width:k+m)][l] + height_flt[k][(l+m>width?l+m-width:l+m)] + height_flt[k][(l-m<0?l-m+width:l-m)])/4;
          height_flt[k][l] = temp + temp2;
        }
      }
    }
    sq_width /= 2;
    max_rand *= REDUCE;
  }
  float tem_max = 0;
  for (unsigned int i = 0; i < height ; i++)
  {
    for (unsigned int j = 0; j < width ; j++)
    {
      tem_max = max(tem_max,fabsf(height_flt[i][j]));
    }
  }
  int adjust = (int)(127.0 / tem_max);
  for (unsigned int i = 0; i < height ; i++)
  {
    for (unsigned int j = 0; j < width ; j++)
    {
      height_map[i][j] = (unsigned char)(tem_max*height_flt[i][j]+128);
    }
    delete height_flt[i];
  }
  return 1;
}

void float_ds_color(float scale = .7)
{
  RGB_heights = new unsigned char**[height];
  for (unsigned int i = 0; i < height; i++)
  {
    RGB_heights[i] = new unsigned char*[width];
    for (unsigned int j = 0; j < width; j++) RGB_heights[i][j] = new unsigned char[3];
  }
  for (int i=0;i<3;i++)
  {
    float_ds(scale);
    for (unsigned int j=0;j<height;j++)
    {
      for (unsigned int k=0;k<width;k++)
      {
        RGB_heights[j][k][i] = height_map[j][k];
      }
    }
  }
  for (unsigned int i = 0; i < height; i++)
  {
    for (unsigned int j = 0; j < width; j++)
    {
      for (int k=0;k<3;k++) outfile << (char)RGB_heights[i][j][k];
      delete RGB_heights[i][j];
    }
    for (unsigned int j = 0; j < buffer_bt; j++) outfile << (char)0;
    delete RGB_heights[i];
  }
}

unsigned long int how_many;

bool water_bodies(int seeds = 1)
{
  Done_Yet = new bool*[height];
  RGB_heights = new unsigned char**[height];
  for (unsigned int i = 0; i < height; i++)
  {
    RGB_heights[i] = new unsigned char*[width];
    Done_Yet[i] = new bool[width];
    for (unsigned int j = 0; j < width; j++) {RGB_heights[i][j] = new unsigned char[3];Done_Yet[i][j]=false;}
  }
  if (!float_ds()) return 0;
  for (unsigned int j=0;j<height;j++)
  {
    for (unsigned int k=0;k<width;k++)
    {
      RGB_heights[j][k][2] = height_map[j][k];
    }
  }
  int ypos,xpos;
  unsigned char depth;
  while (seeds>0)
  {
    ypos = rand()%height;
    xpos = rand()%width;
    depth = height_map[ypos][xpos]+1;
    if ((depth>128)||(depth<96)||(Done_Yet[ypos][xpos])) continue;
    water_fill(xpos,ypos,depth,1);
    if ((float)how_many<pow((float)width,1.25)) continue;
    seeds--;
  }
  for (unsigned int i = 0; i < height; i++)
  {
    for (unsigned int j = 0; j < width; j++)
    {
      if (Done_Yet[i][j]) outfile << (char)(RGB_heights[i][j][0])<<(char)0<<(char)0;
      else outfile << (char)RGB_heights[i][j][2] << (char)RGB_heights[i][j][2] << (char)RGB_heights[i][j][2];
      delete RGB_heights[i][j];
    }
    for (unsigned int j = 0; j < buffer_bt; j++) outfile << (char)0;
    delete RGB_heights[i];
  }
  return 1;
}

void water_fill(int xs, int ys, unsigned char goal, unsigned long int counter = 1)
{
  how_many = max(counter,how_many);
  if (Done_Yet[ys][xs]) return;
  RGB_heights[ys][xs][0]=(255-RGB_heights[ys][xs][2]);
  Done_Yet[ys][xs] = true;
  if ((ys>0)&&((abs(RGB_heights[ys-1][xs][2]-RGB_heights[ys][xs][2])<100)&&(!Done_Yet[ys-1][xs]))&&(RGB_heights[ys-1][xs][2]<goal)) water_fill(xs,ys-1,goal,counter+1);
  if ((ys<(width-1))&&((abs(RGB_heights[ys+1][xs][2]-RGB_heights[ys][xs][2])<100)&&(!Done_Yet[ys+1][xs]))&&(RGB_heights[ys+1][xs][2]<goal)) water_fill(xs,ys+1,goal,counter+1);
  if ((xs>0)&&((abs(RGB_heights[ys][xs-1][2]-RGB_heights[ys][xs][2])<100)&&(!Done_Yet[ys][xs-1]))&&(RGB_heights[ys][xs-1][2]<goal)) water_fill(xs-1,ys,goal,counter+1);
  if ((xs<(width-1))&&((abs(RGB_heights[ys][xs+1][2]-RGB_heights[ys][xs][2])<100)&&(!Done_Yet[ys][xs+1]))&&(RGB_heights[ys][xs+1][2]<goal)) water_fill(xs+1,ys,goal,counter+1);
  return;
}

void basic_avg(float scale = 0.7)
{
  RGB_heights = new unsigned char**[height];
  for (unsigned int i = 0; i < height; i++)
  {
    RGB_heights[i] = new unsigned char*[width];
    for (unsigned int j = 0; j < width; j++) RGB_heights[i][j] = new unsigned char[3];
  }
  for (int i=0;i<3;i++)
  {
    float_ds(scale);
    for (unsigned int j=0;j<height;j++)
    {
      for (unsigned int k=0;k<width;k++)
      {
        RGB_heights[j][k][i] = height_map[j][k];
      }
    }
  }
  for (unsigned int i = 0; i < height; i++)
  {
    for (unsigned int j = 0; j < width; j++)
    {
      outfile << (char)((2*RGB_heights[i][j][0] + RGB_heights[i][j][1])/3);
      outfile << (char)((RGB_heights[i][j][1] + 2*RGB_heights[i][j][2])/3);
      outfile << (char)((RGB_heights[i][j][0] + RGB_heights[i][j][1] + RGB_heights[i][j][2])/3);
      delete RGB_heights[i][j];
    }
    for (unsigned int j = 0; j < buffer_bt; j++) outfile << (char)0;
    delete RGB_heights[i];
  }
}

void adv_avg(float scale = 0.7)
{
  char dummy;
  cout<<"Using floats or characters? (enter c or f) ";
  cin >>dummy;
  int blu1,blu2,blu3,gre1,gre2,gre3,red1,red2,red3,tot1,tot2,tot3;
  cout<<"Enter weight for blue, green, and red:"<<endl;
  cin>>blu1>>gre1>>red1;
  cout<<"Enter weight for blue, green, and red again:"<<endl;
  cin>>blu2>>gre2>>red2;
  cout<<"Enter weight for blue, green, and red one more time:"<<endl;
  cin>>blu3>>gre3>>red3;
  tot1 = blu1+gre1+red1;
  tot2 = blu2+gre2+red2;
  tot3 = blu3+gre3+red3;
  RGB_heights = new unsigned char**[height];
  for (unsigned int i = 0; i < height; i++)
  {
    RGB_heights[i] = new unsigned char*[width];
    for (unsigned int j = 0; j < width; j++) RGB_heights[i][j] = new unsigned char[3];
  }
  for (int i=0;i<3;i++)
  {
    if (dummy=='f') float_ds(scale);
    else if (dummy=='c') diamond_square();
    else {cout<<"Fail! Needs to be f or c."<<endl;return;}
    for (unsigned int j=0;j<height;j++)
    {
      for (unsigned int k=0;k<width;k++)
      {
        RGB_heights[j][k][i] = height_map[j][k];
      }
    }
  }
  for (unsigned int i = 0; i < height; i++)
  {
    for (unsigned int j = 0; j < width; j++)
    {
      outfile << (char)((blu1*RGB_heights[i][j][0] + gre1*RGB_heights[i][j][1] + red1*RGB_heights[i][j][2])/max(tot1,1));
      outfile << (char)((blu2*RGB_heights[i][j][0] + gre2*RGB_heights[i][j][1] + red2*RGB_heights[i][j][2])/max(tot2,1));
      outfile << (char)((blu3*RGB_heights[i][j][0] + gre3*RGB_heights[i][j][1] + red3*RGB_heights[i][j][2])/max(tot3,1));
      delete RGB_heights[i][j];
    }
    for (unsigned int j = 0; j < buffer_bt; j++) outfile << (char)0;
    delete RGB_heights[i];
  }
}

void hght_fill(float scale = 0.7)
{
  RGB_heights = new unsigned char**[height];
  for (unsigned int i = 0; i < height; i++)
  {
    RGB_heights[i] = new unsigned char*[width];
    for (unsigned int j = 0; j < width; j++) RGB_heights[i][j] = new unsigned char[3];
  }
  char dummy;
  cout<<"Using floats or characters? (enter c or f) ";
  cin >>dummy;
  for (int i=0;i<3;i++)
  {
    if (dummy=='f') float_ds(scale);
    else if (dummy=='c') diamond_square();
    else {cout<<"Fail! Needs to be f or c."<<endl;return;}
    for (unsigned int j=0;j<height;j++)
    {
      for (unsigned int k=0;k<width;k++)
      {
        RGB_heights[j][k][i] = height_map[j][k];
      }
    }
  }
  if (dummy=='c')
  {
    unsigned char mini,maxi;
    unsigned char mi[3],ma[3];
    mini = 255;
    maxi = 0;
    for (int k=0;k<3;k++)
    {
      for (int i=0;i<height;i++)
      {
        for (int j=0;j<width;j++)
        {
          mini = min(mini,RGB_heights[i][j][k]);
          maxi = max(maxi,RGB_heights[i][j][k]);
        }
      }
      ma[k]==maxi;
      mi[k]==mini;
    }
    int targ = 0;
    int mark = ma[0]-mi[0];
    if (mark<(ma[1]-mi[1])) {targ = 1; mark = (ma[1]-mi[1]);}
    if (mark<(ma[2]-mi[2])) {targ = 2; mark = (ma[2]-mi[2]);}
    for (unsigned int i = 0; i < height; i++)
    {
      for (unsigned int j = 0; j < width; j++)
      {
        height_map[i][j] = RGB_heights[i][j][targ];
        delete RGB_heights[i][j];
      }
      delete RGB_heights[i];
    }
  }
  else{
  for (unsigned int i = 0; i < height; i++)
  {
    for (unsigned int j = 0; j < width; j++)
    {
      height_map[i][j] = (unsigned char)((RGB_heights[i][j][0] + RGB_heights[i][j][1] + RGB_heights[i][j][2])/3);
      delete RGB_heights[i][j];
    }
    delete RGB_heights[i];
  }}
  unsigned char val;
  for (unsigned int i = 0; i < height; i++)
  {
    for (unsigned int j = 0; j < width; j++)
    {
if(dummy=='f'){
      if (height_map[i][j]>191)
      {
        val = (255 - (height_map[i][j]));
        outfile << (char)val << (char)((val + height_map[i][j])/2) << (char)val;
      }
      else if (height_map[i][j]>127)
      {
        val = (255 - (height_map[i][j]));
        outfile << (char)(val*2) << (char)(val + height_map[i][j]/4) << (char)val;
      }
      else if (height_map[i][j]>63)
      {
        val = (255 - (height_map[i][j] * 2));
        outfile << (char)val << (char)((val + height_map[i][j]/2)) << (char)((val + height_map[i][j]));
      }
      else
      {
        val = (255 - (height_map[i][j]*3));
        outfile << (char)val << (char)val << (char)val;
      }}
else{
      if (height_map[i][j]>130)
      {
        val = (255 - (height_map[i][j]));
        outfile << (char)(val/2) << (char)(val/2) << (char)(val/2);
      }
      else if (height_map[i][j]>110)
      {
        val = (255 - (height_map[i][j]));
        outfile << (char)(val/3) << (char)((val - height_map[i][j]/3)) << (char)(val/3);
      }
      else if (height_map[i][j]>90)
      {
        val = (255 - (height_map[i][j] * 2));
        outfile << (char)val << (char)((val + height_map[i][j]/4)) << (char)((val + height_map[i][j]/2));
      }
      else
      {
        val = (255 - (height_map[i][j]*3));
        outfile << (char)(val/2 + height_map[i][j]/4) << (char)(val/4 + height_map[i][j]/8) << (char)(val>>4);
      }}
    }
    for (unsigned int j = 0; j < buffer_bt; j++) outfile << (char)0;
  }
}

void seeded_dsa()
{
  char dummy;
  int seeds;
  unsigned int wide;
  unsigned int i,j;
  cout<<"Float or short? (Enter f or s) ";
  cin>>dummy;
  cout<<"How many seeds? (number between 0 and 127) ";
  cin>>seeds;
  if (dummy>'r') dummy='s';
  else if (dummy<'g') dummy = 'f';
  else dummy = 's';
  if (seeds<1) seeds = 0;
  else if (seeds>126) seeds = 127;
  unsigned int place = 1;
  unsigned count_it = 0;
  while (place<width) {if (width==(place+1)) break;place = place<<1; count_it++;}
  while (place<height) {if (height==(place+1)) break;place = place<<1; count_it++;}
  wide = place + 1;
  short** highs;
  bool** seeded;
  highs = new short*[wide];
  seeded = new bool*[wide];
  for (i=0;i<wide;i++) {highs[i] = new short[wide];seeded[i]=new bool[wide];for (j=0;j<wide;j++) {highs[i][j] = 0;seeded[i][j]=false;}}
  unsigned int xpos,ypos;
  if (dummy=='s')
  {
    for (i=0;i<seeds;i++)
    {
      xpos = rand()%(wide/2) + (wide/4);
      ypos = rand()%(wide/2) + (wide/4);
      if ((xpos&((wide-2)>>7)) || (ypos&((wide-2)>>7)) || (seeded[ypos][xpos])) {i--; continue;}
      else {highs[ypos][xpos] = ((rand()%8192) - 1024);seeded[ypos][xpos] = true;}
    }
    unsigned int sq_width = wide - 1;
    unsigned int max_rand = 16384;
    for (i=0;i<wide;i++)
    {
      highs[0][i] = -5500; seeded[0][i] = true;
      highs[wide-1][i] = -5500; seeded[wide-1][i] = true;
      highs[i][wide-1] = -5500; seeded[i][wide-1] = true;
      highs[i][0] = -5500; seeded[i][0] = true;
    }
    int x_pos = 0;
    int y_pos = 0;
    short temp;
    short temp2;
    while (sq_width>0)
    {
      for (i = 0; i < (wide - 1); i += sq_width)
      {
        for (j = 0; j < (wide - 1); j += sq_width)
        {
          x_pos = (i + (sq_width/2));
          y_pos = (j + (sq_width/2));
          if (seeded[y_pos][x_pos]);// highs[y_pos][x_pos] = (highs[y_pos][x_pos] + highs[j][i] + highs[j+sq_width][i] + highs[j][i+sq_width] + highs[j+sq_width][i+sq_width])/5 + ((rand()%max_rand - (max_rand/2)));
          else highs[y_pos][x_pos] = (highs[j][i] + highs[j+sq_width][i] + highs[j][i+sq_width] + highs[j+sq_width][i+sq_width])/4 + ((rand()%max_rand - (max_rand/2)));
        }
      }
      int k,l,m;
      m = sq_width/2;
      for (i = 0; i < (wide - 1); i += sq_width)
      {
        for (j = 0; j < (wide - 1); j += sq_width)
        {
          x_pos = (i + m);
          y_pos = (j + m);
        //do up first
          l = x_pos;
          k = j;
          temp2 = rand()%max_rand - max_rand/2;
          if (seeded[k][l]) temp = highs[k][l];// (highs[(k-m<0?k-m+width:k-m)][l] + highs[(k+m>width?k+m-width:k+m)][l] + highs[k][(l+m>width?l+m-width:l+m)] + highs[k][(l-m<0?l-m+width:l-m)] + highs[k][l])/5 + temp2;
          else temp = (highs[(k-m<0?k-m+width:k-m)][l] + highs[(k+m>width?k+m-width:k+m)][l] + highs[k][(l+m>width?l+m-width:l+m)] + highs[k][(l-m<0?l-m+width:l-m)])/4 + temp2;
          highs[k][l] = temp;
        //do left next
          l = i;
          k = y_pos;
          temp2 = rand()%max_rand - max_rand/2;
          if (seeded[k][l]) temp = highs[k][l];// (highs[(k-m<0?k-m+width:k-m)][l] + highs[(k+m>width?k+m-width:k+m)][l] + highs[k][(l+m>width?l+m-width:l+m)] + highs[k][(l-m<0?l-m+width:l-m)] + highs[k][l])/5 + temp2;
          else temp = (highs[(k-m<0?k-m+width:k-m)][l] + highs[(k+m>width?k+m-width:k+m)][l] + highs[k][(l+m>width?l+m-width:l+m)] + highs[k][(l-m<0?l-m+width:l-m)])/4 + temp2;
          highs[k][l] = temp;
        //check and do right then down
          if ((x_pos + sq_width)>width)
          {
            l = i+sq_width;
            k = y_pos;
            highs[k][l] = highs[k][0];
          }
          if ((y_pos + sq_width)>width)
          {
            l = x_pos;
            k = j + sq_width;
            highs[k][l] = highs[0][l];
          }
        }
      }
      sq_width /= 2;
      max_rand = (max_rand>4?max_rand/2:4);
    }
    short placehold = 0;
    short placehol0 = 0;
    short placehol2 = 1;
    short placehol3 = 1;
    for (i=0;i<height;i++)
    {
      for (j=0;j<width;j++)
      {
        placehold = max(highs[i][j],placehold);
        placehol0 = min(highs[i][j],placehol0);
      }
    }
    while ((placehold/placehol2)>127) placehol2++;
    while ((placehol0/placehol3)<(-127)) placehol3++;
    for (i=0;i<height;i++)
    {
      for (j=0;j<width;j++)
      {
        height_map[i][j] = (unsigned char)((highs[i][j]/((highs[i][j]<0)?placehol3:placehol2))+128);
      }
      delete highs[i];
    }
    delete highs;
  }
  else
  {
    float factor;
    cout<<"What scale factor? ";
    cin>>factor;
    if (factor>1) factor = 1;
    else if (factor<0) factor = 0;
    unsigned int sq_width = wide - 1;
    float max_rand = 32;
    float REDUCE = pow(2.0,-factor);
    float** height_flt;
    height_flt = new float*[wide];
    for (i = 0; i < wide; i++)
    {
      height_flt[i] = new float[wide];
      for (j=0;j<wide;j++)
      {
        height_flt[i][j] = 0;
      }
    }
    for (i=0;i<seeds;i++)
    {
      xpos = rand()%wide;
      ypos = rand()%wide;
      if ((xpos&((wide-2)>>9)) || (ypos&((wide-2)>>9)) || (seeded[ypos][xpos]) || (xpos==0) || (ypos==0)) {i--; continue;}
      else {height_flt[ypos][xpos] = (frand()*16 - 4);seeded[ypos][xpos] = true;}
    }
    for (i=0;i<wide;i++)
    {
      height_flt[0][i] = -7;
      height_flt[i][0] = -7;
      height_flt[wide-1][i] = -7;
      height_flt[i][wide-1] = -7;
    }
    height_flt[0][0] = height_flt[0][wide-1] = height_flt[wide-1][0] = height_flt[wide-1][wide-1] = 0;
    int x_pos = 0;
    int y_pos = 0;
    float temp;
    float temp2;
    while(sq_width>0)
    {
      for (i = 0; i < (wide - 1); i += sq_width)
      {
        for (j = 0; j < (wide - 1); j += sq_width)
        {
          x_pos = (i + (sq_width/2));
          y_pos = (j + (sq_width/2));
          if (seeded[y_pos][x_pos]) height_flt[y_pos][x_pos] = height_flt[y_pos][x_pos];// (height_flt[y_pos][x_pos] + height_flt[j][i] + height_flt[j+sq_width][i] + height_flt[j][i+sq_width] + height_flt[j+sq_width][i+sq_width])/5 + frand()*max_rand - max_rand/2;
          else height_flt[y_pos][x_pos] = (height_flt[j][i] + height_flt[j+sq_width][i] + height_flt[j][i+sq_width] + height_flt[j+sq_width][i+sq_width])/4 + frand()*max_rand - max_rand/2;
        }
      }
      int k,l,m;
      m = sq_width/2;
      for (unsigned int i = 0; i < (wide - 1); i += sq_width)
      {
        for (unsigned int j = 0; j < (wide - 1); j += sq_width)
        {
          x_pos = (i + m);
          y_pos = (j + m);
        //do up first
          l = x_pos;
          k = j;
          temp2 = frand()*max_rand - max_rand/2;
          if (seeded[k][l])temp = height_flt[k][l];// (height_flt[(k-m<0?k-m+width:k-m)][l] + height_flt[(k+m>width?k+m-width:k+m)][l] + height_flt[k][(l+m>width?l+m-width:l+m)] + height_flt[k][(l-m<0?l-m+width:l-m)] + height_flt[k][l])/5;
          else temp = (height_flt[(k-m<0?k-m+width:k-m)][l] + height_flt[(k+m>width?k+m-width:k+m)][l] + height_flt[k][(l+m>width?l+m-width:l+m)] + height_flt[k][(l-m<0?l-m+width:l-m)])/4 + temp2;
          height_flt[k][l] = temp;
        //do left next
          l = i;
          k = y_pos;
          temp2 = frand()*max_rand - max_rand/2;
          if (seeded[k][l]) temp = height_flt[k][l];// (height_flt[(k-m<0?k-m+width:k-m)][l] + height_flt[(k+m>width?k+m-width:k+m)][l] + height_flt[k][(l+m>width?l+m-width:l+m)] + height_flt[k][(l-m<0?l-m+width:l-m)] + height_flt[k][l])/5;
          else temp = (height_flt[(k-m<0?k-m+width:k-m)][l] + height_flt[(k+m>width?k+m-width:k+m)][l] + height_flt[k][(l+m>width?l+m-width:l+m)] + height_flt[k][(l-m<0?l-m+width:l-m)])/4 + temp2;
          height_flt[k][l] = temp;
        //check and do right then down
          if ((x_pos + sq_width)>width)
          {
            l = i+sq_width;
            k = y_pos;
            height_flt[k][l] = height_flt[k][0];
          }
          if ((y_pos + sq_width)>width)
          {
            l = x_pos;
            k = j + sq_width;
            height_flt[k][l] = height_flt[0][l];
          }
        }
      }
      sq_width /= 2;
      max_rand *= REDUCE;
    }
    float tem_max = 0;
    float tem_min = 0;
    for (i = 0; i < height ; i++)
    {
      for (j = 0; j < width ; j++)
      {
        tem_max = max(tem_max,height_flt[i][j]);
        tem_min = min(tem_min,height_flt[i][j]);
      }
    }
    float adjust1 = (127.0 / tem_max);
    float adjust2 = (-127.0 / tem_min);
    for (i = 0; i < height ; i++)
    {
      for (j = 0; j < width ; j++)
      {
        height_map[i][j] = (unsigned char)((height_flt[i][j]*(height_flt[i][j]<0?adjust2:adjust1))+128);
      }
      delete height_flt[i];
    }
    delete height_flt;
  }
  unsigned char val;
  for (i = 0; i < height; i++)
  {
    for (j = 0; j < width; j++)
    {
      if (height_map[i][j]>130)
      {
        val = (255-(height_map[i][j]));
        outfile << (char)(val/2) << (char)(val/2) << (char)(val/2);
      }
      else if (height_map[i][j]>110)
      {
        val = (255 - (height_map[i][j]));
        outfile << (char)(val/3) << (char)((val - height_map[i][j]/3)) << (char)(val/3);
      }
      else if (height_map[i][j]>90)
      {
        val = (255 - (height_map[i][j] * 2));
        outfile << (char)val << (char)((val + height_map[i][j]/4)) << (char)((val + height_map[i][j]/2));
      }
      else
      {
        val = (255 - (height_map[i][j]*5)/2);
        outfile << (char)(140-(val/2 + height_map[i][j]/4)) << (char)(78-(val/4 + height_map[i][j]/8)) << (char)(val>>4);
      }
    }
    for (unsigned int j = 0; j < buffer_bt; j++) outfile << (char)0;
  }
}

void TwoD_MD()
{
  float horizontal[width];
  float vertical[height];
  unsigned int i,j;
  for (i=0;i<height;i++) vertical[i]=0;
  for (i=0;i<width;i++)  horizontal[i]=0;
  unsigned int seg_width = width-1;
  float max_rand = 32.0;
  unsigned int xpos;
  while (seg_width>0)
  {
    for (i=0;i<width-1;i+=seg_width)
    {
      xpos = i+seg_width/2;
      horizontal[xpos] = (frand()*max_rand - max_rand/2.0) + (horizontal[i]+horizontal[i+seg_width])/2.0;
    }
    max_rand /= 2.0;
    seg_width /= 2;
  }
  seg_width = height-1;
  max_rand = 32.0;
  while (seg_width>0)
  {
    for (i=0;i<height-1;i+=seg_width)
    {
      xpos = i+seg_width/2;
      vertical[xpos] = (frand()*max_rand - max_rand/2.0) + (vertical[i]+vertical[i+seg_width])/2.0;
    }
    max_rand /= 2.0;
    seg_width /= 2;
  }
  float tem_max,tem_min;
  tem_max=tem_min=0;
  for (i = 0;i<width;i++)
  {
    for (j = 0;j<height; j++)
    {
      tem_max = max(tem_max,horizontal[i]*vertical[j]);
      tem_min = min(tem_min,horizontal[i]*vertical[j]);
    }
  }
  float adjust1,adjust2;
  adjust1 = 127.0 / tem_max;
  adjust2 = -127.0 / tem_min;
  for (i=0;i<height;i++)
  {
    for (j=0;j<width;j++)
    {
      height_map[i][j] = (unsigned char)((horizontal[j]*vertical[i]*(horizontal[j]*vertical[i]<=0?adjust2:adjust1))+128);
    }
  }
  unsigned char val;
  for (i = 0; i < height; i++)
  {
    for (j = 0; j < width; j++)
    {
      if (height_map[i][j]>230)
      {
        outfile << (char)(height_map[i][j]) << (char)(height_map[i][j]) << (char)(height_map[i][j]);
      }
      else if (height_map[i][j]>184)
      {
        val = (255-(height_map[i][j]));
        outfile << (char)(val) << (char)(val) << (char)(val);
      }
      else if (height_map[i][j]>156)
      {
        val = (255 - height_map[i][j]+50);
        outfile << (char)(val/4) << (char)((val - (height_map[i][j]/3))) << (char)(val/4);
      }
      else if (height_map[i][j]>128)
      {
        val = (255 - (height_map[i][j] * 7)/5);
        outfile << (char)val << (char)((val + height_map[i][j]/4)) << (char)((val + height_map[i][j]/2));
      }
      else
      {
        val = (255 - (height_map[i][j]*2));
        outfile << (char)(140-(val/3 + height_map[i][j]/5)) << (char)(78-(val/5 + height_map[i][j]/8)) << (char)(val>>4);
      }
   }
    for (unsigned int j = 0; j < buffer_bt; j++) outfile << (char)0;
  }
}
