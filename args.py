class Args:
    def __init__(self, output, contigs, plastic, mappings):
        # Check if all arguments are strings and not empty or contains spaces
        invalid_args = [name for name, value in locals().items() if not (isinstance(value, str) and value and ' ' not in value) and name != 'self']
        if invalid_args:
            raise ValueError(f"Please provide valid arguments for: {', '.join(invalid_args)}")
        
        self.output = output
        self.contigs = contigs
        self.plastic = plastic.lower()
        self.mappings = mappings

