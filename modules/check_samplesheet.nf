@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

// List of allowed extensions: h5ad, RDS
allowed_extensions = ['h5ad', 'rds']

def check_samplesheet(input_file) {
    def samples = []
    for(line in parseCsv(new FileReader(input_file))) {
        adata_path = line['input_adata']
        adata = file(adata_path)
        assert adata.exists() : "File does not exist: ${adata_path}"
        extension = adata_path.split('\\.').last().toLowerCase()
        // Make sure extension is allowed
        assert extension in allowed_extensions : "File extension not allowed: ${adata_path}"
        meta = line.toMap()
        meta.remove('input_adata')
        samples << [meta, adata, extension]
    }
    return samples
}
