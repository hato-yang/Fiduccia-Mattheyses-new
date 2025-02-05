struct array_metadata
{
	int tolerance;
	int total_area;
	int number_of_cells;
	int number_of_nets;
	struct net **NET_array;
	// For some reason, the struct fails to free unless this dummy_variable is not assigned a value
	//   I believe this is true because it follows the struct net**, although I don't know why
	int dummy_variable;
	struct cell **CELL_array;
};

struct are_metadata
{
	long total_area;
	long total_area2;
	long die_area;
	int tolerance;
	int tolerance1; // 儲存第二製程最大面積
	int tolerance_Macro;
	int tolerance1_Macro; // 儲存第二製程最大面積
	// 新增 u 數據的欄位
	int u_value1;
	int u_value2;
	long u_value3;
};

struct condensed *read_in_data_to_arrays(char *, char *);

int count_cells_in_are_file(char *);
struct are_metadata *read_in_are_file(struct cell **, char *);

int count_nets_in_netD_file(char *);
int read_in_netD_file(struct cell **, struct net **, char *);
