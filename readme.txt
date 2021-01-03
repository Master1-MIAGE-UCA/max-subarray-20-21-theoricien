parallel:
	void montee (EMODE, struct tablo *, struct tablo *);
	void descente (struct tablo *,struct tablo *);
	void final (struct tablo *, struct tablo *);
	void max_montee (struct tablo *, struct tablo *);
	void max_descente (struct tablo *, struct tablo *);
	void max_final (struct tablo *,struct tablo *);
	void find_max_subarray (struct tablo *, struct tablo *);
	int main (int, char *[]);

non parallel:
	void _debug (const char *, ...);
	void sum (EMODE, struct tablo *, struct tablo **);
	void generate_array (struct tablo **);
	void print_array (struct tablo *);
	void read_array (struct tablo **, char *);
	void max (struct tablo *, struct tablo **);	
	size_t new_size (int);
	size_t file_size (FILE *);
	struct tablo * init_tablo (size_t);
