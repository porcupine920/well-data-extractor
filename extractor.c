#include <math.h>
#include <stdio.h>

#define PI 3.14159265358979

static int data_num;
static int gt_len;
static float r_jd = 0.02f;
static float r_num;
static float reactive[3000]; //reactive
static float current[3000]; //current
static float voltage[3000]; //voltage
static float power[3000]; //active
static float displacement[3000];
static float load[3000];
static float factor[3000];
static float data[3000];

void get_dft(float a[], int len, float pcos[], float psin[], int order)
{
	int i, j;
	float t = 2.0f * PI / len;
	float c = cos(t);
	float s = sin(t);
	float s1 = 0, c1 = 1, u0;
	t = 2.0f / len;
	for (i = 0; i < order; ++i) {
		float u1 = 0.0f;
		float u2 = 0.0f;
		for (j = len - 1; j > 0; --j) {
			u0 = a[j] + 2 * c1 * u1 - u2;
			u2 = u1;
			u1 = u0;
		}
		pcos[i] = t * (a[0] + u1 * c1 - u2);
		psin[i] = t * u1 * s1;
		u0 = c * c1 - s * s1;
		s1 = c * s1 + s * c1;
		c1 = u0;
	}
}

void tran_data()
{
	int i, j, n;
	float i1[gt_len], u1[gt_len], p1[gt_len], s1[gt_len], g1[gt_len];
	float x, ljys, min;
	float psin[7], pcos[7];
	get_dft(displacement, gt_len, pcos, psin, 6);
	min = 10;
	for (i = -80; i < 80; ++i) {
		x = i * 2 * PI / gt_len;
		ljys = 0;
		for (j = 1; i < 6; ++j)
			ljys += 2 * j * (psin[j] * cos(j * x) - pcos[j] * sin(j * x));
		if (fabs(ljys) < min) {
			min = fabs(ljys);
			n = i;
		}
	}
	for (i = 0; i < gt_len; ++i) {
		i1[i] = current[i];
		u1[i] = voltage[i];
		p1[i] = power[i];
		s1[i] = displacement[i];
		g1[i] = factor[i];
	}
	for (i = 0; i < gt_len; ++i) {
		j = (i + n + gt_len) % gt_len;
		current[i] = i1[j];
		voltage[i] = u1[j];
		power[i] = p1[j];
		displacement[i] = s1[j];
		factor[i] = g1[j];
	}
}

void put_dft(float a[], int num, float pcos[], float psin[], int order)
{
	int i, j;
	float t = 2 * PI / num;
	for (i = 0; i < num; ++i) {
		a[i] = pcos[0] / 2;
		for (j = 1; j < order; ++j)
			a[i] += pcos[j] * cos(i * j * t) + psin[j] * sin(i * j * t);
	}
}

void decode_wave(float data[], int len, float data2[])
{
	int avgdot, i, d_num;
	float q_sum;

	float data1[len];
	
	avgdot = 20;
	q_sum = 0;
	d_num = 0;
	for (i = 0; i < avgdot; ++i) {
		q_sum = q_sum + data[i];
		d_num++;
		data1[i] = q_sum / d_num;
	}

	for (i = avgdot; i < len; ++i) {
		q_sum -= data[i - avgdot] + data[i];
		data1[i] = q_sum / d_num;
	}
	len -= avgdot;
	for (i = 0; i < len; ++i)
		data2[i] = data1[i + avgdot];
}

void get_cycle(char id[], float data[], int len)
{
	int t1, t2, t3, tn, found_num, t_num;
	int dat_flg[3000];
	float data2[3000];

	decode_wave(data, len, data2);
	int q_max = -999, i, j;
	for (i = 0; i < len; ++i)
		if (data2[i] > q_max)
			q_max = data2[i];
	float q_num = 0.96f;
	do {
		float q_max1 = q_max * q_num;
		found_num = 0;
		for (i = 0; i < len; ++i) {
			if (data2[i] < q_max1)
				dat_flg[i] = 0;
			else {
				dat_flg[i] = 1;
				found_num++;
			}
		}
		if (found_num / len < 0.08) {
			q_num -= 0.02;
			if (q_num < 0.5) {
				gt_len = len;
				return;
			}
		}
	} while (found_num / len < 0.08);

	for (i = 0; i < len; ++i) {
		if (dat_flg[i] == 0)
			continue;
		t1 = i;
		dat_flg[i] = 0;
		for (j = i + 1; j < len; ++j) {
			if (dat_flg[j] == 0) {
				t2 = j;
				j = i = len;
			} else {
				dat_flg[j] = 0;
			}
		}
	}

	t1 = (int)((t1 + t2 - 1) / 2);
	dat_flg[t1] = 1;
	int last_slope = 99, lock_flg = 0, last_t = 0;
	for (t3 = t2 + 190; t3 < t2 + (len - t2 + 1) / 2; ++t3) {
		if (dat_flg[t3]) {
			t_num = t3 - t1;
			tn = (int)((len - t1) / t_num - 1);
			float x1sum = 0, y1sum = 0, x2sum = 0, y2sum = 0, xysum = 0;
			for (j = t1; j < t1 + t_num; ++j) {
				x1sum += data2[j];
				y1sum += data2[j + tn * t_num];
				x2sum += data2[j] * data2[j];
				y2sum += data2[j + tn * t_num] * data2[j + tn * t_num];
				xysum += data2[j] * data2[j + tn * t_num];
			}
			float lxx = x2sum - powf(x1sum, 2.0f) / t_num;
			float lyy = y2sum - powf(y1sum, 2.0f) / t_num;
			float lxy = xysum - x1sum * y1sum / t_num;

			float b_iu = lxy / lxx;
			float a_iu = (y1sum - b_iu * x1sum) / t_num;
			float r_iu = lxy / sqrtf(fabsf(lxx * lyy));

			double slope = fabs(1.0 - r_iu);
			if (slope <= r_jd) {
				last_slope = slope;
				last_t = t_num;
			} else {
				t_num = last_t + 1;
				gt_len = t_num;
				r_num = r_iu;
				return;
			}
		}
	}
	t_num = len;
}

void tx_org(int flag)
{
	float psin[2], pcos[2];
	float mydata[data_num];
	if (flag == 1) {
		//generate file path and load data
	}
	r_num = 0;
	gt_len = 3000;
	get_cycle("first_displacement", displacement, data_num);

	get_dft(displacement, gt_len, pcos, psin, 1);
	int i, t_num;
	if (gt_len > 2000 || gt_len <= 250 || sqrt(pow(psin[1], 2) + pow(pcos[1], 2)) < 1) {
		r_num = 0;
		gt_len = 3000;
		get_cycle("first_displacement", power, data_num);
		if (gt_len > 2000 || gt_len <= 250) {
			for (i = 0; i < data_num - 500; ++i)
				mydata[i] = power[500 + i];
			r_num = 0;
			get_cycle("second_power", mydata, (data_num - 500));
			if (gt_len > 2000) {
				gt_len = data_num;
				t_num = data_num;
			}
		}
	} else
	    tran_data();
	int cycle_dot = gt_len;
	int iu_num = gt_len;
	/*
	 * UI related operation and function invocation
	 */
	 float uu, ii, pp2, ii2;
	 float p3 = 0, ppl = 0;
	 int pi;
	 for (i = 0; i < gt_len; ++i) {
		 pi = (i + 1) % gt_len;
		 ppl = (factor[pi] + factor[i]) / 2 * (displacement[pi] - displacement[i]);
		 p3 += ppl;
	 }
	 // UI related operation

	 float pp_sum = 0;
	 for (i = 0; i < gt_len; ++i)
		 pp_sum += power[i];
	 // set value of UI component
	 float i1_max = 0;
	 for (i = 0; i < gt_len / 2; ++i)
		 if (current[i] > i1_max)
			 i1_max = current[i];
	 // set value of UI component
	 float i2_max = 0;
	 for (i = gt_len / 2; i < gt_len; ++i)
		 if (current[i] > i2_max)
			 i2_max = current[i];
	 // set value of UI component
	 i1_max = (int)(i1_max * 10 + 0.5f) / 10;
	 i2_max = (int)(i2_max * 10 + 0.5f) / 10;
	 // set value of UI component

	 float p1_avg = 0;
	 for (i = 0; i < gt_len / 2; ++i)
		 p1_avg += power[i];
	 int np1_avg = (int)((p1_avg / i * 10.0f + 0.5f) / 10);

	 float p2_avg = 0;
	 for (; i < gt_len; ++i)
		 p2_avg += power[i];
	 int np2_avg = (int)((p2_avg / (gt_len / 2) * 10 + 0.5f) / 10);
	 // display the quotient of np1_avg and np2_avg

	 // averge voltage, current, apparent power
	 float u_avg = 0, i_avg = 0, ap = 0, i_sq = 0, p_sq = 0;
	 for (i = 0; i < gt_len; ++i) {
		 u_avg += voltage[i];
		 i_avg += current[i];
		 ap += voltage[i] * current[i];
		 i_sq += current[i] * current[i];
		 p_sq += power[i] * power[i];
	 }
	 // set values on UI accordingly

}

int main(int argc, char* argv[])
{
	FILE *f = fopen(argv[1], "r");
	fclose(f);
	return 0;
}
