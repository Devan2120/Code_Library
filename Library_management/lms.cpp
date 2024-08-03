#include<unistd.h>
#include<bits/stdc++.h>

using namespace std;

void clrscr(){                               //clear screen
    system("cls");                           //windows
    //system("clear");                         //linux
}

void exi(int cond) {   
    if(cond==0){
        cout<<"\nExiting...";
        sleep(1);
        exit(0);
    }
}

void file_error(){
    cout<<"\n\n\tFile Opening Error";
    exi(0);
}

void prints(char str[]){                                 //print string
    for(int p=0;p<strlen(str);p++){
        if(p==30){
            cout<<"-\n";
        }
        cout<<str[p];
    }
}

void ret(){                                             
    cout<<"\n\nPress Any Button";
    cin.ignore();
}

void langc(){                                      //language choice
    cout<<"\nENTER THE LANGUAGE OF THE BOOK:";
    cout<<(" _________________");
    cout<<"\n\n";
    cout<<("\t\t\tLANGUAGE CHOICE");
    cout<<"\n\n\n\t\t";
    cout<<("1.ENGLISH");
    cout<<"\n\t\t";
    cout<<("2.MALAYALAM") ;
    cout<<"\n\t\t";
    cout<<("3.HINDI");
    cout<<"\n\t\t";
    cout<<("4.OTHERS");
    cout<<"\n\n\t\t";
    cout<<("PLEASE ENTER THE NUMBER TO CHOOSE THE CORRESPONDING OPTION or ENTER 0 TO QUIT: ");
}

void genrec(){
    cout<<"\nENTER THE GENRE OF BOOK";
    cout<<(" _________________");
    cout<<"\n\n";
    cout<<("\t\t\tGENRE CHOICE");
    cout<<"\n\n\n\t\t";
    cout<<("1.MAGAZINES/COMICS");
    cout<<"\n\t\t";
    cout<<("2.NOVEL");
    cout<<"\n\t\t";
    cout<<("3.CHEMISTRY");
    cout<<"\n\t\t";
    cout<<("4.PHYSICS");
    cout<<"\n\t\t";
    cout<<("5.MATHS");
    cout<<"\n\t\t";
    cout<<("6.HISTORY");
    cout<<"\n\t\t";
    cout<<("7.LANGUAGE");
    cout<<"\n\t\t";
    cout<<("8.BIOLOGY");
    cout<<"\n\t\t";
    cout<<("9.GEOGRAPHY");
    cout<<"\n\t\t";
    cout<<("10.COMPUTER");
    cout<<"\n\t\t";
    cout<<("11.COMMERCE");
    cout<<"\n\t\t";
    cout<<("12.OTHERS");
    cout<<("\n\n\t\tPLEASE ENTER THE NUMBER TO CHOOSE THE CORRESPONDING OPTION or ENTER 0 TO QUIT: ");    
}


class book{					 //class book
    int copy;
    char bno[12];
    int lang;
    char bname[40];
    int genre;
    char aname[20];

    public:
    
    void gen_bookno();
    void create_bk(){ 												//enter book details
	    cout<<"\nENTER THE NAME OF THE BOOK: ";
		gets(bname);
		cout<<"\nENTER THE NAME OF THE AUTHOR: ";
		gets(aname);
		langc();
        cin>>lang;
        exi(lang);
        clrscr();
        genrec();
        cin>>genre;
        exi(genre);
        gen_bookno();
        cout<<"\n\t\tTHE BOOK NO IS: ";
        puts(bno);
        cout<<"\n\t\tENTER THE NUMBER OF COPIES: ";
        cin>>copy;
        cin.ignore();
        ret();
	}

    void display_bk(){									//display book details     
	    cout<<"\n\t\tBOOK NUMBER IS :";
		puts(bno);
		cout<<"\n\t\tBOOK NAME IS :";
		puts(bname);
		cout<<"\n\t\tAUTHOR NAME IS:";
		puts(aname);
		cout<<"\n\t\tCOPIES AVAILABLE :";
		cout<<copy;
        //cin.ignore();
	}

    char* get_bno(){
        return bno;
    }
    char* get_bname(){
        return bname;
    }
    int get_copy(){
        return copy;
    }
    void i_copy(){
        copy=copy+1;
    }
    void d_copy(){
        copy=copy-1;
    }

    void report(){                             //details of a book
        cout<<"\n";
        cout<<bno<<"\t";
        prints(bname);
        cout<<"\t";
        if(strlen(bname)<15){
            cout<<"\t\t";
        }
        prints(aname);
        if(strlen(aname)<10){
            cout<<"\t";
        }
        cout<<"\t"<<copy<<"\n";
    }
};

void append_bk(book b){                      //new book to file
    ofstream out;
    out.open("book.dat",ios::app);
    if(!out){
        file_error();
    }
    else{
        out.write((char*)&b,sizeof(b));
        out.close();
    }    
}

string noc(int no){                           //integer to string of 3 charecters
    char str[3];
    itoa(no,str,10);
    int precision = 3 - strlen(str);
    return string(precision, '0').append(str);
}

void book::gen_bookno(){                       //fn to generate book no;
    string bno_temp;
    int id;
    bno_temp=noc(lang)+"/"+noc(genre)+"/";
    fstream inout;
    inout.open("id.dat",ios::in);
    if(!inout){
        inout.close();
        inout.open("id.dat",ios::out);
        id=1;
        inout<<id;
    }
    else{
        inout>>id;
        id++;
        inout.close();
        inout.open("id.dat",ios::trunc|ios::out);
        inout<<id;  
    }
    inout.close();
    bno_temp+=noc(id);
    strcpy(bno,bno_temp.c_str());
}    

void modify_bk(){                 //modify data in file
    book bk,obj;
    fstream inout("book.dat",ios::in|ios::out);
    if(!inout){
        file_error();
    }
    else{
        char no[12];
        cout<<"\n\t\tEnter The Book No: ";
        gets(no);
        int flag=0;
        int count=0;
        while(!inout.eof() && flag==0){
            count++;
            inout.read((char*)&obj,sizeof(obj));
            if(!strcmp(obj.get_bno(),no)){
                cout<<"\n\n\t\t The Current Details Are: \n";
                obj.display_bk();
                cout<<"\n\t\t Enter the New Details: \n";
                bk.create_bk();
                inout.seekp((count-1)*sizeof(obj));
                inout.write((char*)&bk,sizeof(bk));
                flag=1;
                cout<<"\n\tRecord Modified";
            }
        }
        if(flag==0){
                cout<<"\n\tNo Such Record Exit";
            }
        inout.close();
    }        
    ret();
}

void remove_bk(){                    //removes a book
    book bk;
    int flag=0;
    ifstream inout("book.dat",ios::in);
    ofstream inout_temp("temp.dat",ios::out);
    if(!inout || !inout_temp){
        file_error();
    }
    else{
        char no[12];
        cout<<"\n\t\tEnter The Book No:";
        gets(no);
        while(!inout.eof()){
            inout.read((char*)&bk,sizeof(bk));
            if(strcmp(bk.get_bno(),no)){
                inout_temp.write((char*)&bk,sizeof(bk));
            }
            else{
                flag=1;
                cout<<"\n\t The Deleted Record is: ";
                bk.display_bk(); 
            }
        }
        if(flag==0){
            cout<<"\n\t No Such Record Exist";
        }
        inout.close();
        inout_temp.close();
        char filename[]="book.dat";
        remove(filename);
        rename("temp.dat","book.dat");
    }
    ret();    
}

void search_b_name(){                                        //search book by name
    char ename[40];
	cout<<"\n\tEnter The Book name To be Searched: ";
	gets(ename);
	int found=0;
	fstream infile;
	book bk;
	infile.open("book.dat",ios::in);
    if(!infile){
        file_error();
    }
    else{
        while(!infile.eof()&&found==0)
            {	infile.read((char*)&bk,sizeof(bk));
                if(!strcmpi(bk.get_bname(),ename))
                {	bk.display_bk();
                    found=1;
                }
            }
        if(found==0)
        cout<<"\n\tSuch A Book Does Not Exist";
        infile.close();
    } 
    ret();   
}

void displayall_bk(){
    ifstream infile;
	book bk;
	clrscr();
	infile.open("book.dat",ios::in);
	if(!infile){
        file_error();
	}
    else{
        cout<<"\n\n\t\tBook List\n\n";
        cout<<"================================================================================\n";
        cout<<"Book Number"<<setw(15)<<"Book Name"<<setw(20)<<"Author"<<setw(15)<<"Copies\n";
        cout<<"================================================================================\n";
        while(infile.read((char*)&bk,sizeof(book))){
            bk.report();
        }
        infile.close();
    }
    ret();    
}

void search_b_no(){                                          //search book by no 
    char bnum[12];
	cout<<"\n\tEnter The Book no To be Searched: ";
	gets(bnum);
	int found=0;
	fstream infile;
	book bk;
	infile.open("book.dat",ios::in);
    if(!infile){
        file_error();
    }
    else{
        while(!infile.eof()&&found==0)
            {	infile.read((char*)&bk,sizeof(bk));
                if(!strcmpi(bk.get_bno(),bnum))
                {	bk.display_bk();
                    found=1;
                }
            }
        if(found==0)
        cout<<"\n\tBook Does Not Exist";
        infile.close();
    }
    ret();
}

void book_list(){
    int choice;
    do{
        clrscr();
        cout<<("\n\t\t\tBOOK LIST");
        cout<<"\n\n";
        cout<<("\n\t\t1.Search By Name");
	    cout<<("\n\t\t2.Search By Number");
	    cout<<("\n\t\t3.Display All Books");
        cout<<("\n\t\t4.MAIN MENU");
        cout<<("\n\n\tPlease Enter The Number To Choose The Corresponding Option Or Enter 0 To Quit: ");
        cin>>choice;
        cin.ignore();
        switch(choice){
            case 0: {exi(0); break;}

            case 1: {search_b_name(); break;}

            case 2: {search_b_no(); break;}

            case 3: {displayall_bk(); break;}

            case 4: {break;}

            default: {cout<<"\n\t\tInvalid Choice. Try Again Or Enter 0 To Quit";}
        }
        clrscr();
    }while(choice!=4);
}

class student{                          //class for students
    int adno;
    char name[20];
    char std[2];
    char div[2];                              
    char issued_book[12];
    int token;

    public:

    void create_st();
    void display_st();

    int get_adno(){
        return adno;
    }
    int get_token(){
        return token;
    }
    void i_token(){
        token=1;
    }
    void d_token(){
        token=0;
    }
    char* get_book(){
        return issued_book;
    }
    void book_issued(char bk[]){
        strcpy(issued_book,bk);
    }
    char* get_name(){
        return name;
    }
    void report(){                                    //prints student record
        cout<<"\n\n\t" ;
        cout<<adno<<"\t";
        prints(name);
        if(strlen(name)<10){
            cout<<"\t";
        }
        cout<<"\t\t";
        prints(std);
        cout<<"-";
        prints(div);
        cout<<"\t";
        if(token==1){
            cout<<"\tYes";
        }
        else{
            cout<<"\tNo";
        }
        cout<<endl;
   }
};

void append_st(student st){                   //adds data to file
    ofstream out;
	out.open("student.dat",ios::app);
    if(!out){
        file_error();
    }
    else{
		out.write((char*)&st,sizeof(st));
		out.close();
    }    
}

void student::create_st(){                           //creates a new student
    cout<<"\n ENTER ADMISSION NO: ";
	cin>>adno;
    cin.ignore();
	cout<<"\n ENTER THE NAME OF THE STUDENT: ";
	gets(name);
	cout<<"\n ENTER THE STANDARD OF THE STUDENT: ";
	gets(std);
	cout<<"\n ENTER THE DIVISION OF THE STUDENT: ";
	gets(div);
    strcpy(issued_book,"Null");
    token=0;
    ret();
}

void modify_st(){                                   //modify exixting student data
    student st,obj;
    fstream inout("student.dat",ios::in|ios::out);
    if(!inout){
        file_error();
    }
    else{
        cout<<"\n\t\tEnter The Admission No Of The Student: ";
        int no;
        cin>>no;
        int flag=0;
        int count=0;
        while(!inout.eof() && flag==0){
            count++;
            inout.read((char*)&obj,sizeof(obj));
            if(obj.get_adno()==no){
                cout<<"\n\n\t\t The Current Details Are: \n";
                obj.display_st();
                cout<<"\n\t\t Enter the New Details: \n";
                st.create_st();
                st.book_issued(obj.get_book());
                if(obj.get_token()==1){
                    st.i_token();
                }
                inout.seekp((count-1)*sizeof(obj));
                inout.write((char*)&st,sizeof(st));
                flag=1;
                cout<<"\n\tRecord Modified";
            }
        }
        if(flag==0){
                cout<<"\n\tNo Such Record Exit";
            }
        inout.close();
    }        
    ret();
}

void student::display_st(){                                     //display student data
    cout<<"\nADMISSION NUMBER IS: "<<adno;
	cout<<"\nNAME OF THE STUDENT IS: ";
	puts(name);
	cout<<"\nCLASS: ";
	puts(std);
	cout<<"\nDIV: ";
	puts(div);
	cout<<"\nBOOKS ISSUED: ";
    if(token==1){
            cout<<"\tYes";
        }
    else{
        cout<<"\tNo";
    }
    //cin.ignore();
}

void remove_st(){                                  //removes a student
    student st;
    int flag=0;
    ifstream inout("student.dat",ios::in);
    ofstream inout_temp("temp.dat",ios::out);
    if(!inout || !inout_temp){
        file_error();
    }
    else{
        cout<<"\n\t\tEnter The Admission No Of The Student: ";
        int no;
        cin>>no;
        while(!inout.eof()){
            inout.read((char*)&st,sizeof(st));
            if(st.get_adno()!=no){
                inout_temp.write((char*)&st,sizeof(st));
            }
            else{
                flag=1;
                if(st.get_token()==1){
                    cout<<"\n\t Cannot Delete Record: Pending Book";
                    inout_temp.write((char*)&st,sizeof(st));
                }
                else{
                    cout<<"\n\t The Deleted Record is: ";
                    st.display_st();
                } 
            }
        }
        if(flag==0){
            cout<<"\n\t No Such Record Exist";
        }
        inout.close();
        inout_temp.close();
        char filename[]="student.dat";
        remove(filename);
        rename("temp.dat","student.dat");
    } 
    ret();   
}

void search_st_no(){                                   //search student by admission no
    int stnum;
	cout<<"\n\tEnter The Student Admission no To be Searched: ";
	cin>>stnum;
	int found=0;
	fstream infile;
	student st;
	infile.open("student.dat",ios::in);
    if(!infile){
        file_error();
    }
    else{
        while(!infile.eof()&&found==0)
            {	infile.read((char*)&st,sizeof(st));
                if(stnum==st.get_adno())
                {	st.display_st();
                    found=1;
                }
            }
        if(found==0)
        cout<<"\nSTUDENT DOES NOT EXIST";
        infile.close();
    }
    ret();
}

void search_st_name(){                             //search student by name
    char ename[20];
	cout<<"\n\tEnter The Student name To be Searched: ";
	gets(ename);
	int found=0;
	fstream infile;
	student st;
	infile.open("student.dat",ios::in);
    if(!infile){
        file_error();
    }
    else{
        while(!infile.eof()&&found==0)
            {	infile.read((char*)&st,sizeof(st));
                if(!strcmpi(st.get_name(),ename))
                {	st.display_st();
                    found=1;
                }
            }
        if(found==0)
        cout<<"\nSuch A Student Does Not Exist";
        infile.close();
    }
    ret();    
}

void displayall_st(){                              //display all students
    ifstream infile;
	student st;
	clrscr();
	infile.open("student.dat",ios::in);
	if(!infile){
        file_error();
	}
    else{
        cout<<"\n\n\t\tStudent List\n\n";
        cout<<"================================================================================\n";
        cout<<"Admission Number"<<setw(15)<<"Student Name"<<setw(10)<<"Class"<<setw(15)<<"Book Issued\n";
        cout<<"================================================================================\n";
        while(infile.read((char*)&st,sizeof(st))){
            st.report();
        }
        infile.close();
    }    
    ret();
}

void student_list(){
    int choice;
    do{
        clrscr();
        cout<<("\n\t\t\tSTUDENTS LIST");
        cout<<"\n\n";
        cout<<("\n\t\t1.Search By Name");
	    cout<<("\n\t\t2.Search By Number");
	    cout<<("\n\t\t3.Display All Students");
        cout<<("\n\t\t4.MAIN MENU");
        cout<<("\n\n\tPlease Enter The Number To Choose The Corresponding Option Or Enter 0 To Quit: ");
        cin>>choice;
        cin.ignore();
        switch(choice){
            case 0: {exi(0); break;}

            case 1: {search_st_name(); break;}

            case 2: {search_st_no(); break;}

            case 3: {displayall_st(); break;}

            case 4: {break;}

            default: {cout<<"\n\t\tInvalid Choice. Try Again Or Enter 0 To Quit";}
        }
    }while(choice!=4);
}


void pend_book(){ 						//display pending books
    student st;
    book bk;
    int sn,bn;
    fstream fs,fb;
    fs.open("student.dat",ios::in);
    fb.open("book.dat",ios::in);
    if(!fs || !fb){
        file_error();
    }
    else{
        int flag=0;
        cout<<"\n\t\t  PENDING BOOKS: ";
        cout<<"\n\t Student Name"<<"\t Book Name"<<"\t Book Number";
        while(!fs.eof()){
            fs.read((char*)&st,sizeof(st));
            if(st.get_token()==1){
                while(!fb.eof()){
                    fb.read((char*)&bk,sizeof(bk));
                    if(!strcmpi(bk.get_bno(),st.get_book())){
                        flag=1;
                        cout<<"\n\t   "<<st.get_name();
                        if(sizeof(st.get_name())<15){
                            cout<<"\t";
                        }
                        cout<<"\t  "<<bk.get_bname();
                        
                        cout<<"\t "<<st.get_book();
                    }
                }
            }
        }
        if(flag==0){
            cout<<"\n\tNone";
        }    
        fs.close();
        fb.close();
    }
    ret();
}

void book_deposit(){                   //fn to return a book
    student st;
    book bk;
    int count1=0;
    int count2=0;
    int flag1=0;
    int flag2=0;
    char bn[12];
    int an;
    cout<<"\n\t Enter Book no: ";
    gets(bn);
    cout<<"\n\t Enter Student Admission no: ";
    cin>>an;
    fstream fb,fs;
    fs.open("student.dat",ios::in|ios::out);
    fb.open("book.dat",ios::in|ios::out);
    if(!fs || !fb){
        file_error();
    }
    else{
        while(!fs.eof() && flag1==0){
            fs.read((char*)&st,sizeof(st));
            count1++;
            if((an==st.get_adno())&&(!strcmpi(bn,st.get_book()))){
                flag1=1;
                while(!fb.eof() && flag2==0){
                    count2++;
                    fb.read((char*)&bk,sizeof(bk));
                    if(!strcmpi(bn,bk.get_bno())){
                        flag2=1;
                        st.d_token();
                        bk.i_copy();
                        char null[]="Null";
                        st.book_issued(null);
                        fs.seekp((count1-1)*sizeof(st),ios::beg);
                        fs.write((char*)&st,sizeof(st));
                        fb.seekp((count2-1)*sizeof(bk),ios::beg);
                        fb.write((char*)&bk,sizeof(bk));
                        cout<<"\n\tBook Returned";
                        break;
                    }            
                }
            }    
        }
        if(flag1==0 || flag2==0){
            cout<<"\n\tBook Returning Error";
        }
    }
    cin.ignore();
    ret();
}

void book_issue(){
    student st;
    book bk;
    int count1=0;
    int count2=0;
    int flag1=0;
    int flag2=0;
    char bn[12];
    int an;
    cout<<"\n\t Enter Book no:";
    gets(bn);
    cout<<"\n\t Enter Student Admission no: ";
    cin>>an;
    fstream fb,fs;
    fs.open("student.dat",ios::in|ios::out);
    fb.open("book.dat",ios::in|ios::out);
    if(!fs || !fb){
        file_error();
    }
    else{
        while(!fs.eof() && flag1==0){
            fs.read((char*)&st,sizeof(st));
            count1++;
            if((an==st.get_adno())&&(!strcmp("Null",st.get_book()))){
                flag1=1;
                while(!fb.eof() && flag2==0){
                    count2++;
                    fb.read((char*)&bk,sizeof(bk));
                    if(!strcmp(bn,bk.get_bno())){
                        flag2=1;
                        st.i_token();
                        bk.d_copy();
                        st.book_issued(bn);
                        fs.seekp((count1-1)*sizeof(st),ios::beg);
                        fs.write((char*)&st,sizeof(st));
                        fb.seekp((count2-1)*sizeof(bk),ios::beg);
                        fb.write((char*)&bk,sizeof(bk));
                        cout<<"\n\tBook Issued";
                        break;
                    }            
                }
            }    
        }
        if(flag1==0 || flag2==0){
            cout<<"\n\tBook Issuing Error";
        }
    }
    cin.ignore();
    ret();
}

int enter_pass(){
    int epass;
    cout<<"\n\n";
    cout<<("+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|+|");
    cout<<"\n";
    cout<<"\n\tEnter Password\n";
    cin>>epass;
    return epass;
}

int check_pass(){
    clrscr();
    int epass;
    int pass;
    fstream inout;
    inout.open("pass.dat",ios::in|ios::out);
    if(!inout){
        inout.open("pass.dat",ios::out);
        cout<<"\n\tSetting up Your Libary Management Software\n";
        cout<<"\n\tLet's First Set a Password\n";
        epass=enter_pass();
        inout<<epass;
        return 1;
    }
    else{
        cout<<"\n\n\tEnter The Password To Login";
        epass=enter_pass();
        inout>>pass;
        if(pass==epass){
            return 1;
        }
        else{
            return 0;
        }    
    }
    inout.close();
}

void change_pass(){
    cout<<"\n\t\tPasword change";
    cout<<"\n\tEnter Old Password";
    int epass=enter_pass();
    fstream inout;
    inout.open("pass.dat",ios::in);
    if(!inout){
        file_error();
    }
    else{
        int pass;
        inout>>pass;
        if(pass==epass){
            inout.close();
            cout<<"\n\tEnter New Password";
            int n_pass=enter_pass();
            inout.open("pass.dat",ios::trunc|ios::out);
            inout<<n_pass;
            cout<<"\n\tPassword Changed";
        }
        else{
            cout<<"\n\t Wrong Password";
        }
        inout.close();
    }
    cin.ignore();
    ret();
}


int main(){
    book bk;
    student st;
    int choice;
    int flag=check_pass();
    if(flag==1){
        do{
	   		clrscr();
	   		cout<<("________________");
	   		cout<<("\t\tW.E.L.C.O.M.E");
	   		cout<<("\t\t  ________________");
	  		cout<<"\n\n";
			cout<<("\t\t\t\tMAIN MENU");
			cout<<"\n\n\n";
			cout<<"\t";
			cout<<("1.ADD A BOOK");
			cout<<"\t\t";
			cout<<("2.BOOK LIST");
			cout<<"\t\t";
			cout<<("3.MODIFY A BOOK");
			cout<<"\n\n\t";
			cout<<("4.DELETE A BOOK");
			cout<<"\t\t";
			cout<<("5.ADD NEW STUDENT");
			cout<<"\t";
			cout<<("6.DELETE A STUDENT");
			cout<<"\n\n\t";
			cout<<("7.MODIFY A STUDENT");
			cout<<"\t";
			cout<<("8.STUDENTS LIST");
			cout<<"\t\t";
			cout<<("9.BOOK ISSUE");
			cout<<"\n\n\t";
			cout<<("10.BOOK RETURN");
            cout<<"\t\t";
			cout<<"11.PENDING BOOKS";
			cout<<"\t";
			cout<<"12.CHANGE PASSWORD";
			cout<<"\n\n\n";
			cout<<(" PLEASE ENTER THE NUMBER TO CHOOSE THE CORRESPONDING OPTION or ENTER 0 TO QUIT: ");
            cin>>choice;
            cin.ignore();
            clrscr();
            switch(choice){

                case 0: {exi(0); break;}

                case 1: {bk.create_bk(); append_bk(bk); break;}

                case 2: {book_list(); break;}

                case 3: {modify_bk(); break;}

                case 4: {remove_bk(); break;}

                case 5: {st.create_st(); append_st(st); break;}

                case 6: {remove_st(); break;}

                case 7: {modify_st(); break;}

                case 8: {student_list(); break;}

                case 9: {book_issue(); break;}

                case 10: {book_deposit(); break;}

                case 11: {pend_book(); break;}

                case 12: {change_pass(); break;}

                default: {cout<<"\n\t\tInvalid Choice. Try Again Or Enter 0 To Quit";}
            }
        }while(choice!=0);
    }    
    else{
        cout<<"\n\tWrong Password";
        exi(0);
    }
    return 0;
}

















