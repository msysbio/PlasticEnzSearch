import flet as ft
from args import Args
import PlasticEnz
import io
import sys
import traceback
import logging
import csv

class PlasticEnzymeSearch:
    def __init__(self, page: ft.Page):
        # Initialize the page and set its title
        self.page = page
        self.page.title = "Plastic Enzyme Search"

        # Define the types of plastics
        self.plastic_types = ('all', 'pbat', 'nylon', 'ab-hydrolase', 'pet', 'pbsa', 'pha', 'pcl', 'phb', 'pla', 'cutinase')

        # Create text fields to display selected files/directory
        self.selected_contigs = ft.TextField()
        self.selected_mappings = ft.TextField()
        self.directory_path = ft.TextField()

        # Create a text widget to display selected plastics
        self.selected_plastics_text = ft.TextField()

        # Create markdown widget to display the output
        self.output_markdown = ft.Markdown()

        # Create progress ring widget
        self.progress_ring = ft.Container(content=ft.ProgressRing(), alignment=ft.alignment.center)

        # Create submit button and add it to the page
        self.submit_btn = ft.ElevatedButton(text="Run PlastEnzSearch", on_click=self.button_clicked)

        # Center the widgets
        self.centered_submit_btn = ft.Container(content=self.submit_btn, alignment=ft.alignment.center)
        self.centered_output_markdown = ft.Container(content=self.output_markdown, alignment=ft.alignment.center)

    # Define result handlers for file pickers
    def pick_contigs_result(self, e: ft.FilePickerResultEvent):
        # Update the selected contigs text field with the selected files
        self.selected_contigs.value = ", ".join(map(lambda f: f.path, e.files)) if e.files else "No files selected!"
        self.selected_contigs.update()

    def pick_mappings_result(self, e: ft.FilePickerResultEvent):
        # Update the selected mappings text field with the selected files
        self.selected_mappings.value = ", ".join(map(lambda f: f.path, e.files)) if e.files else "No files selected!"
        self.selected_mappings.update()

    def get_directory_result(self, e: ft.FilePickerResultEvent):
        # Update the directory path text field with the selected directory
        self.directory_path.value = e.path if e.path else "No directory selected!"
        self.directory_path.update()

    # Define the click handler for plastic types
    def check_item_clicked(self, e, plastic):
        # Toggle the checkbox value
        e.control = not e.control

        # If 'all' is selected, uncheck all other checkboxes
        if plastic == 'all':
            for checkbox in self.plastic_checkboxes:
                checkbox.value = False
            e.control = True

        # Update the selected plastics text field with the selected plastics
        selected_plastics = [item.label for item in self.plastic_checkboxes if item.value]
        if 'all' in selected_plastics and len(selected_plastics) > 1:
            selected_plastics = ['all']
        self.selected_plastics_text.value = ", ".join(selected_plastics)
        self.selected_plastics_text.update()
        self.page.update()

    # Define click handler for the submit button
    def button_clicked(self, e):
        # Redirect standard output to a string buffer
        old_stdout = sys.stdout
        sys.stdout = buffer = io.StringIO()

        try:
            # Create an Args object with the selected values
            args = Args(
                output=self.directory_path.value,
                contigs=self.selected_contigs.value,
                plastic=self.selected_plastics_text.value,
                mappings=self.selected_mappings.value
            )

            # Add progress ring to the page
            self.page.add(self.progress_ring)

            # Run the main function
            PlasticEnz.main(args)

            # Remove progress ring from the page
            self.page.remove(self.progress_ring)

            # Read the TSV file
            plastics = []
            if args.plastic == 'all':
                plastics = self.plastic_types
            else:
                plastics = args.plastic.split(',')
            for plastic in plastics:

                try:
                    with open(f'{args.output}/temps/{plastic}/mapping_summary.tsv', 'r') as f:
                        reader = csv.reader(f, delimiter='\t')
                        data = list(reader)
                    # Create a DataTable and add it to the page
                    table = ft.DataTable(items=data)
                    self.page.append(table)
                except FileNotFoundError:
                    pass

        except Exception as e:
            # Write the error message to the buffer
            buffer.write(str(e))
            # Log the full traceback
            logging.error(traceback.format_exc())

        # Reset standard output
        sys.stdout = old_stdout

        # Update markdown widget with the output
        self.output_markdown.value = buffer.getvalue()
        self.output_markdown.update()
        self.page.update()

    def run(self):
        # Create file pickers
        pick_contigs_dialog = ft.FilePicker(on_result=self.pick_contigs_result)
        pick_mappings_dialog = ft.FilePicker(on_result=self.pick_mappings_result)
        get_directory_dialog = ft.FilePicker(on_result=self.get_directory_result)

        # Add file pickers to the page overlay
        self.page.overlay.extend([pick_contigs_dialog, pick_mappings_dialog, get_directory_dialog])

        # Create checkboxes for each plastic type
        self.plastic_checkboxes = [
            ft.Checkbox(label=plastic, on_change=lambda e, plastic=plastic: self.check_item_clicked(e, plastic))
            for plastic in self.plastic_types
        ]

        # Add on_change event to each checkbox
        for checkbox in self.plastic_checkboxes:
            checkbox.on_change = lambda e: self.check_item_clicked(e, checkbox.label)

        # Create a row for plastic types
        plastic_row = ft.Row(self.plastic_checkboxes, spacing=8)

        # Add file pickers, plastic types row, and selected plastics text to the page
        self.page.add(
            ft.Row([ft.Text("Output Directory:"), self.directory_path]),
            ft.Row([
                ft.ElevatedButton(
                    "Choose",
                    icon=ft.icons.FOLDER_OPEN,
                    on_click=lambda _: get_directory_dialog.get_directory_path(),
                ),
            ]),
            ft.Row([ft.Text("Contigs Files:"), self.selected_contigs]),
            ft.Row([
                ft.ElevatedButton(
                    "Pick Files",
                    icon=ft.icons.UPLOAD_FILE,
                    on_click=lambda _: pick_contigs_dialog.pick_files(allow_multiple=True),
                ),
            ]),
            ft.Row([ft.Text("Mappings Files:"), self.selected_mappings]),
            ft.Row([
                ft.ElevatedButton(
                    "Pick Files",
                    icon=ft.icons.UPLOAD_FILE,
                    on_click=lambda _: pick_mappings_dialog.pick_files(allow_multiple=True),
                ),
            ]),
            ft.Row([ft.Text("Plastic types:"), self.selected_plastics_text]),
            plastic_row,
            ft.Divider(height=10, color="white"),
            self.centered_submit_btn, 
            self.centered_output_markdown
        )

if __name__ == '__main__':
    # Define a function to run the app
    def run_app(page: ft.Page):
        # Create an instance of PlasticEnzymeSearch with the page argument
        PlasticEnzymeSearch(page).run()
    # Pass run_app as the target function to ft.app
    # ft.app will call run_app and pass the page argument to it
    ft.app(target=run_app)#, view=ft.WEB_BROWSER)
