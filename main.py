import flet as ft
from args import Args
import PlasticEnz
import io
import sys
import traceback
import logging
import os
import webbrowser

class PlasticEnzymeSearch:
    """A class representing the PlasticEnz GUI.

    This class provides functionality to run the Plastic Tool with a graphical user interface.

    Args:
        page (ft.Page): The page object representing the user interface.

    Attributes:
        page (ft.Page): The page object representing the user interface.
        plastic_types (tuple): A tuple containing the types of plastics.
        selected_contigs (ft.TextField): A text field to display selected contig files.
        selected_mappings (ft.TextField): A text field to display selected mapping files.
        directory_path (ft.TextField): A text field to display the selected directory path.
        selected_plastics_text (ft.TextField): A text field to display the selected plastic types.
        output_markdown (ft.Markdown): A markdown widget to display the output.
        progress_ring (ft.Container): A container widget containing a progress ring.
        submit_btn (ft.ElevatedButton): A button to submit the search.
        centered_submit_btn (ft.Container): A container widget to center the submit button.
        centered_output_markdown (ft.Container): A container widget to center the output markdown.

    """

    def __init__(self, page: ft.Page):
        """Initialize the PlasticEnzymeSearch object.

        Args:
            page (ft.Page): The page object representing the user interface.

        """
        # Initialize the page and set its title
        self.page = page
        self.page.title = "Plastic Enzyme Search"

        # Define the types of plastics
        self.plastic_types = ('all', 'pbat', 'nylon', 'ab-hydrolase', 'pet', 'pbsa', 'pha', 'pcl', 'phb', 'pla', 'cutinase')

        textfield_width = 1000
        # Create text fields to display selected files/directory
        self.selected_contigs = ft.TextField(width=textfield_width)
        self.selected_mappings = ft.TextField(width=textfield_width)
        self.directory_path = ft.TextField(width=textfield_width)

        # Create a text widget to display selected plastics
        self.selected_plastics_text = ft.TextField(width=textfield_width)

        # Create markdown widget to display the output
        self.output_markdown = ft.Markdown()

        # Create progress ring widget
        self.progress_ring = ft.Container(content=ft.ProgressRing(), alignment=ft.alignment.center)

        # Create submit button and add it to the page
        self.submit_btn = ft.ElevatedButton(text="Run PlastEnzSearch", on_click=self.button_clicked)

        # Center the widgets
        self.centered_submit_btn = ft.Container(content=self.submit_btn, alignment=ft.alignment.center)
        self.centered_output_markdown = ft.Container(content=self.output_markdown, alignment=ft.alignment.center)

    def pick_contigs_result(self, e: ft.FilePickerResultEvent):
        """Handle the result of the contigs file picker.

        This method updates the selected contigs text field with the selected files.

        Args:
            e (ft.FilePickerResultEvent): The event object containing the selected files.

        """
        # Update the selected contigs text field with the selected files
        self.selected_contigs.value = ",".join(map(lambda f: f.path, e.files)) if e.files else "No files selected!"
        self.selected_contigs.update()

    def pick_mappings_result(self, e: ft.FilePickerResultEvent):
        """Handle the result of the mappings file picker.

        This method updates the selected mappings text field with the selected files.

        Args:
            e (ft.FilePickerResultEvent): The event object containing the selected files.

        """
        # Update the selected mappings text field with the selected files
        self.selected_mappings.value = ",".join(map(lambda f: f.path, e.files)) if e.files else "No files selected!"
        self.selected_mappings.update()

    def get_directory_result(self, e: ft.FilePickerResultEvent):
        """Handle the result of the directory picker.

        This method updates the directory path text field with the selected directory.

        Args:
            e (ft.FilePickerResultEvent): The event object containing the selected directory.

        """
        # Update the directory path text field with the selected directory
        self.directory_path.value = e.path if e.path else "No directory selected!"
        self.directory_path.update()

    def check_item_clicked(self, e, plastic):
        """Handle the click event of a plastic type checkbox.

        This method toggles the checkbox value and updates the selected plastics text field.

        Args:
            e: The event object containing the checkbox information.
            plastic (str): The plastic type associated with the checkbox.

        """
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
        self.selected_plastics_text.value = ",".join(selected_plastics)
        self.selected_plastics_text.update()
        self.page.update()




    def button_clicked(self, e):
        """Handle the click event of the submit button.

        This method runs the PlasticEnzSearch application and updates the output markdown widget.

        Args:
            e: The event object containing the button information.

        """
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

            # Create the html report button and add it to the page
            self.launch_url_btn = ft.ElevatedButton(text="html report", on_click=lambda e, args=args: webbrowser.open(args.output+'/abundances.html'))
            self.page.add(ft.Container(content=self.launch_url_btn, alignment=ft.alignment.center))

        except Exception as e:
            # Write the error message to the buffer
            buffer.write(str(e))
            # Log the full traceback
            logging.error(traceback.format_exc())

        # Reset standard output
        sys.stdout = old_stdout

        # Update markdown widget with the output
        #self.output_markdown.value = buffer.getvalue()
        #self.output_markdown.update()
        self.page.update()

    def auto_button_clicked(self, e: ft.FilePickerResultEvent):
        """Handle the click event of the auto button.

        This method automatically fills in the text fields based on the selected directory.

        Args:
            e (ft.FilePickerResultEvent): The event object containing the selected directory.

        """

        # Check if the directory path is valid
        if e.path:
            # Create the output directory path
            output_directory_path = os.path.join(e.path, "output")
            self.directory_path.value = output_directory_path
            self.directory_path.update()

            # Create the output directory if it does not exist
            if not os.path.exists(output_directory_path):
                os.makedirs(output_directory_path)

            # Get all files in the directory that end with ".fa" and update the selected contigs text field
            contigs_files = [os.path.join(e.path, file) for file in os.listdir(e.path) if file.endswith(".fa")]
            self.selected_contigs.value = ",".join(contigs_files) if contigs_files else "No files selected!"
            self.selected_contigs.update()

            # Get all files in the directory that end with ".bam" and update the selected mappings text field
            mappings_files = [os.path.join(e.path, file) for file in os.listdir(e.path) if file.endswith(".bam")]
            self.selected_mappings.value = ",".join(mappings_files) if mappings_files else "No files selected!"
            self.selected_mappings.update()

            # Update the selected plastics text field with "all"
            self.selected_plastics_text.value = "all"
            self.selected_plastics_text.update()


    def run(self):
        """Run the Plastic Enzyme Search application.

        This method sets up the user interface and runs the application.

        """
        # Create file pickers
        auto_dialog = ft.FilePicker(on_result=self.auto_button_clicked)
        pick_contigs_dialog = ft.FilePicker(on_result=self.pick_contigs_result)
        pick_mappings_dialog = ft.FilePicker(on_result=self.pick_mappings_result)
        get_directory_dialog = ft.FilePicker(on_result=self.get_directory_result)

        # Add file pickers to the page overlay
        self.page.overlay.extend([auto_dialog, pick_contigs_dialog, pick_mappings_dialog, get_directory_dialog])

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
            ft.Row([
                ft.ElevatedButton(
                    "Auto Fill",
                    icon=ft.icons.FOLDER_OPEN,
                    on_click=lambda _: auto_dialog.get_directory_path(),
                ),
            ]),
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
    def run_app(page: ft.Page):
        """Run the Plastic Enzyme Search application.

        This function creates an instance of the PlasticEnzymeSearch class and runs the application.

        Args:
            page (ft.Page): The page object representing the user interface.

        """
        PlasticEnzymeSearch(page).run()

    ft.app(target=run_app)
    #ft.app(target=run_app, view=ft.AppView.WEB_BROWSER)
